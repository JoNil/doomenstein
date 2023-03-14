use glam::{vec2, Mat2, Vec2};
use minifb::{Key, Scale, ScaleMode, Window, WindowOptions};
use parse_display::FromStr;
use std::{
    f32::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU},
    fs,
    time::Duration,
};

const SCREEN_WIDTH: usize = 384;
const SCREEN_HEIGHT: usize = 216;

const EYE_Z: f32 = 1.65;
const HFOV: f32 = 120.0 * std::f32::consts::PI / 180.0;
const FRAC_HFOV_2: f32 = HFOV / 2.0;
const VFOV: f32 = 0.5;

const ZNEAR: f32 = 0.0001;
const ZFAR: f32 = 128.0;

fn ifnan(x: f32, alt: f32) -> f32 {
    if x.is_nan() {
        alt
    } else {
        x
    }
}

// see: https://en.wikipedia.org/wiki/Line–line_intersection
// compute intersection of two line segments, returns None if there is
// no intersection.
fn intersect_segs(a0: Vec2, a1: Vec2, b0: Vec2, b1: Vec2) -> Option<Vec2> {
    let a_diff = a0 - a1;
    let b_diff = b0 - b1;
    let d = a_diff.perp_dot(b_diff);

    if d.abs() < 0.000001 {
        return None;
    }

    let ab = a0 - b0;
    let t = ab.perp_dot(b_diff) / d;
    let u = ab.perp_dot(a_diff) / d;
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        Some(a0 - t * a_diff)
    } else {
        None
    }
}

fn abgr_mul(col: u32, a: u32) -> u32 {
    let br = ((col & 0xFF00FF) * a) >> 8;
    let g = ((col & 0x00FF00) * a) >> 8;
    0xFF000000 | (br & 0xFF00FF) | (g & 0x00FF00)
}

#[derive(Copy, Clone, Default, FromStr)]
#[display("{a.x} {a.y} {b.x} {b.y} {portal}")]
struct Wall {
    #[from_str(default)]
    a: Vec2,
    #[from_str(default)]
    b: Vec2,
    portal: usize,
}

// sector id for "no sector"
const SECTOR_NONE: usize = 0;
const SECTOR_MAX: usize = 128;

#[derive(Copy, Clone, Default, FromStr)]
#[display("{_id} {firstwall} {nwalls} {zfloor} {zceil}")]
struct Sector {
    _id: i32,
    firstwall: usize,
    nwalls: usize,
    zfloor: f32,
    zceil: f32,
}

#[derive(Default)]
struct Level {
    sectors: [Sector; 32],
    sectors_count: usize,

    walls: [Wall; 32],
    walls_count: usize,
}

#[derive(Copy, Clone, Default)]
struct Camera {
    pos: Vec2,
    angle: f32,
    camera_rot: Mat2,
    sector: usize,
}

fn world_to_camera(camera: &Camera, p: Vec2) -> Vec2 {
    camera.camera_rot.mul_vec2(p - camera.pos)
}

// convert angle in [-(HFOV / 2)..+(HFOV / 2)] to X coordinate
fn screen_angle_to_x(angle: f32) -> usize {
    ((SCREEN_WIDTH as f32 / 2.0)
        * (1.0 - (((angle + (FRAC_HFOV_2)) / HFOV) * FRAC_PI_2 - FRAC_PI_4).tan())) as usize
}

// noramlize angle to +/-PI
fn normalize_angle(a: f32) -> f32 {
    a - (TAU * ((a + PI) / TAU).floor())
}

// load sectors from file -> state
fn load_level(path: &str) -> Result<Level, String> {
    // sector 0 does not exist
    let mut level = Level {
        sectors_count: 1,
        ..Default::default()
    };

    let str = fs::read_to_string(path).map_err(|e| format!("Unable to read file {path:?}: {e}"))?;

    enum State {
        Sector,
        Wall,
    }

    let mut parse_state = State::Sector;

    for (i, line) in str.lines().enumerate() {
        let line = line.trim();

        if line == "[SECTOR]" {
            parse_state = State::Sector;
            continue;
        }
        if line == "[WALL]" {
            parse_state = State::Wall;
            continue;
        }

        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        match parse_state {
            State::Sector => {
                let sector = line
                    .parse()
                    .map_err(|e| format!("Unable to parse sector {path}:{}: {e}", i + 1))?;
                level.sectors[level.sectors_count] = sector;
                level.sectors_count += 1;
            }
            State::Wall => {
                let wall = line
                    .parse()
                    .map_err(|e| format!("Unable to parse wall {path}:{}: {e}", i + 1))?;
                level.walls[level.walls_count] = wall;
                level.walls_count += 1;
            }
        }
    }

    Ok(level)
}

fn verline(pixels: &mut [u32], x: usize, y0: usize, y1: usize, abgr: u32) {
    for y in y0..=y1 {
        pixels[(SCREEN_HEIGHT - y - 1) * SCREEN_WIDTH + x] = abgr.to_be().rotate_right(8);
    }
}

// -1 right, 0 on, 1 left
fn point_side(p: Vec2, a: Vec2, b: Vec2) -> f32 {
    -(p - a).perp_dot(b - a)
}

// point is in sector if it is on the left side of all walls
fn point_in_sector(level: &Level, sector: &Sector, p: Vec2) -> bool {
    for i in 0..sector.nwalls {
        let wall = &level.walls[sector.firstwall + i];

        if point_side(p, wall.a, wall.b) > 0.0 {
            return false;
        }
    }

    true
}

fn render(pixels: &mut [u32], level: &Level, camera: &Camera) {
    let mut y_lo = [0; SCREEN_WIDTH];
    let mut y_hi = [SCREEN_HEIGHT - 1; SCREEN_WIDTH];

    // track if sector has already been drawn
    let mut sectdraw = [false; SECTOR_MAX];

    // calculate edges of near/far planes (looking down +Y axis)
    let zdl = Mat2::from_angle(FRAC_HFOV_2).mul_vec2(Vec2::Y);
    let zdr = Mat2::from_angle(-FRAC_HFOV_2).mul_vec2(Vec2::Y);
    let znl = zdl * ZNEAR;
    let znr = zdr * ZNEAR;
    let zfl = zdl * ZFAR;
    let zfr = zdr * ZFAR;

    let mut queue = vec![(camera.sector, 0, SCREEN_WIDTH - 1)];

    while let Some((id, x0, x1)) = queue.pop() {
        if sectdraw[id] {
            continue;
        }

        sectdraw[id] = true;

        let sector = &level.sectors[id];

        for i in 0..sector.nwalls {
            let wall = &level.walls[sector.firstwall + i];

            // translate relative to player and rotate points around player's view
            let op0 = world_to_camera(camera, wall.a);
            let op1 = world_to_camera(camera, wall.b);

            // wall clipped pos
            let (mut cp0, mut cp1) = (op0, op1);

            // both are negative -> wall is entirely behind player
            if cp0.y <= 0.0 && cp1.y <= 0.0 {
                continue;
            }

            // angle-clip against view frustum
            let mut ap0 = normalize_angle(f32::atan2(cp0.y, cp0.x) - FRAC_PI_2);
            let mut ap1 = normalize_angle(f32::atan2(cp1.y, cp1.x) - FRAC_PI_2);

            // clip against view frustum if both angles are not clearly within
            // HFOV
            if cp0.y < ZNEAR || cp1.y < ZNEAR || ap0 > FRAC_HFOV_2 || ap1 < -FRAC_HFOV_2 {
                let il = intersect_segs(cp0, cp1, znl, zfl);
                let ir = intersect_segs(cp0, cp1, znr, zfr);

                // recompute angles if points change
                if let Some(il) = il {
                    cp0 = il;
                    ap0 = normalize_angle(f32::atan2(cp0.y, cp0.x) - FRAC_PI_2);
                }

                if let Some(ir) = ir {
                    cp1 = ir;
                    ap1 = normalize_angle(f32::atan2(cp1.y, cp1.x) - FRAC_PI_2);
                }
            }

            if ap0 < ap1 {
                continue;
            }

            if (ap0 < -FRAC_HFOV_2 && ap1 < -FRAC_HFOV_2)
                || (ap0 > FRAC_HFOV_2 && ap1 > FRAC_HFOV_2)
            {
                continue;
            }

            // "true" xs before portal clamping
            let tx0 = screen_angle_to_x(ap0);
            let tx1 = screen_angle_to_x(ap1);

            // bounds check against portal window
            if tx0 > x1 {
                continue;
            }
            if tx1 < x0 {
                continue;
            }

            let wall_diff = wall.b - wall.a;
            let wallshade = 16.0 * f32::sin(f32::atan2(wall_diff.x, wall_diff.y)) + 1.0;

            let x0 = usize::clamp(tx0, x0, x1);
            let x1 = usize::clamp(tx1, x0, x1);

            let z_floor = sector.zfloor;
            let z_ceil = sector.zceil;
            let nz_floor = if wall.portal > 0 {
                level.sectors[wall.portal].zfloor
            } else {
                0.0
            };
            let nz_ceil = if wall.portal > 0 {
                level.sectors[wall.portal].zceil
            } else {
                0.0
            };

            let sy0 = ifnan((VFOV * SCREEN_HEIGHT as f32) / cp0.y, 1e10);
            let sy1 = ifnan((VFOV * SCREEN_HEIGHT as f32) / cp1.y, 1e10);

            let yf0 = (SCREEN_HEIGHT as i32 / 2) + ((z_floor - EYE_Z) * sy0) as i32;
            let yc0 = (SCREEN_HEIGHT as i32 / 2) + ((z_ceil - EYE_Z) * sy0) as i32;
            let yf1 = (SCREEN_HEIGHT as i32 / 2) + ((z_floor - EYE_Z) * sy1) as i32;
            let yc1 = (SCREEN_HEIGHT as i32 / 2) + ((z_ceil - EYE_Z) * sy1) as i32;
            let nyf0 = (SCREEN_HEIGHT as i32 / 2) + ((nz_floor - EYE_Z) * sy0) as i32;
            let nyc0 = (SCREEN_HEIGHT as i32 / 2) + ((nz_ceil - EYE_Z) * sy0) as i32;
            let nyf1 = (SCREEN_HEIGHT as i32 / 2) + ((nz_floor - EYE_Z) * sy1) as i32;
            let nyc1 = (SCREEN_HEIGHT as i32 / 2) + ((nz_ceil - EYE_Z) * sy1) as i32;
            let txd = tx1 - tx0;
            let yfd = yf1 - yf0;
            let ycd = yc1 - yc0;
            let nyfd = nyf1 - nyf0;
            let nycd = nyc1 - nyc0;

            for x in x0..=x1 {
                let shade = if x == x0 || x == x1 {
                    192
                } else {
                    255 - wallshade as u32
                };

                // calculate progress along x-axis via tx{0,1} so that walls
                // which are partially cut off due to portal edges still have
                // proper heights
                let xp = ifnan((x - tx0) as f32 / txd as f32, 0.0);

                // get y coordinates for this x
                let tyf = (xp * yfd as f32) as i32 + yf0;
                let tyc = (xp * ycd as f32) as i32 + yc0;
                let yf = i32::clamp(tyf, y_lo[x] as i32, y_hi[x] as i32) as usize;
                let yc = i32::clamp(tyc, y_lo[x] as i32, y_hi[x] as i32) as usize;

                // floor
                if yf > y_lo[x] {
                    verline(pixels, x, y_lo[x], yf, 0xFFFF0000);
                }

                // ceiling
                if yc < y_hi[x] {
                    verline(pixels, x, yc, y_hi[x], 0xFF00FFFF);
                }

                if wall.portal > 0 {
                    let tnyf = (xp * nyfd as f32) as i32 + nyf0;
                    let tnyc = (xp * nycd as f32) as i32 + nyc0;
                    let nyf = i32::clamp(tnyf, y_lo[x] as i32, y_hi[x] as i32) as usize;
                    let nyc = i32::clamp(tnyc, y_lo[x] as i32, y_hi[x] as i32) as usize;

                    verline(pixels, x, nyc, yc, abgr_mul(0xFF00FF00, shade));
                    verline(pixels, x, yf, nyf, abgr_mul(0xFF0000FF, shade));

                    y_hi[x] = usize::clamp(yc.min(nyc).min(y_hi[x]), 0, SCREEN_HEIGHT - 1);
                    y_lo[x] = usize::clamp(yf.max(nyf).max(y_lo[x]), 0, SCREEN_HEIGHT - 1);
                } else {
                    verline(pixels, x, yf, yc, abgr_mul(0xFFD0D0D0, shade));
                }
            }

            if wall.portal > 0 {
                queue.push((wall.portal, x0, x1));
            }
        }
    }
}

fn main() {
    let mut window = Window::new(
        "raycast-rs",
        SCREEN_WIDTH,
        SCREEN_HEIGHT,
        WindowOptions {
            scale: Scale::X4,
            scale_mode: ScaleMode::Stretch,
            ..Default::default()
        },
    )
    .unwrap();

    window.limit_update_rate(Some(Duration::from_micros(16667)));

    let mut pixels = vec![0; SCREEN_WIDTH * SCREEN_HEIGHT * 4];

    let level = load_level("level.txt").unwrap();

    let mut camera = Camera {
        pos: vec2(3.0, 3.0),
        angle: 0.0,
        sector: 1,
        ..Default::default()
    };

    while window.is_open() && !window.is_key_down(Key::Escape) {
        let rot_speed = 3.0 * 0.016;
        let move_speed = 3.0 * 0.016;

        let mut new_camera = camera;

        if window.is_key_down(Key::Right) {
            new_camera.angle -= rot_speed;
        }

        if window.is_key_down(Key::Left) {
            new_camera.angle += rot_speed;
        }

        new_camera.camera_rot = Mat2::from_angle(-new_camera.angle);

        if window.is_key_down(Key::Up) {
            new_camera.pos += new_camera.camera_rot.row(1) * move_speed;
        }

        if window.is_key_down(Key::Down) {
            new_camera.pos += new_camera.camera_rot.row(1) * -move_speed;
        }

        // update player sector
        {
            // BFS neighbors in a circular queue, player is likely to be in one
            // of the neighboring sectors

            let mut queue = vec![new_camera.sector];

            {
                let sector = &level.sectors[new_camera.sector];

                for j in 0..sector.nwalls {
                    let wall = &level.walls[sector.firstwall + j];

                    if wall.portal > 0 {
                        queue.push(wall.portal)
                    }
                }
            }

            let mut found = SECTOR_NONE;

            while let Some(id) = queue.pop() {
                let sector = &level.sectors[id];

                if point_in_sector(&level, sector, new_camera.pos) {
                    found = id;
                    break;
                }
            }

            if found == SECTOR_NONE {
                println!("Player Collided with wall!");
            } else {
                new_camera.sector = found;
                camera = new_camera;
            }
        }

        render(&mut pixels, &level, &camera);
        window
            .update_with_buffer(&pixels, SCREEN_WIDTH, SCREEN_HEIGHT)
            .unwrap();
    }
}
