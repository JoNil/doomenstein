use minifb::{Key, Scale, ScaleMode, Window, WindowOptions};
use parse_display::FromStr;
use std::{
    collections::VecDeque,
    f32::{
        consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU},
        NAN,
    },
    fs::{self},
    time::Duration,
};

const SCREEN_WIDTH: usize = 384;
const SCREEN_HEIGHT: usize = 216;

const EYE_Z: f32 = 1.65;
const HFOV: f32 = 90.0 * std::f32::consts::PI / 180.0;
const FRAC_HFOV_2: f32 = HFOV / 2.0;
const VFOV: f32 = 0.5;

const ZNEAR: f32 = 0.0001;
const ZFAR: f32 = 128.0;

#[derive(Copy, Clone, Default)]
struct V2 {
    x: f32,
    y: f32,
}

fn v2(x: f32, y: f32) -> V2 {
    V2 { x, y }
}

fn ifnan(x: f32, alt: f32) -> f32 {
    if x.is_nan() {
        alt
    } else {
        x
    }
}

// -1 right, 0 on, 1 left
fn point_side(p: V2, a: V2, b: V2) -> f32 {
    -(((p.x - a.x) * (b.y - a.y)) - ((p.y - a.y) * (b.x - a.x)))
}

// rotate vector v by angle a
fn rotate(v: V2, a: f32) -> V2 {
    v2(
        (v.x * a.cos()) - (v.y * a.sin()),
        (v.x * a.sin()) + (v.y * a.cos()),
    )
}

// see: https://en.wikipedia.org/wiki/Line–line_intersection
// compute intersection of two line segments, returns (NAN, NAN) if there is
// no intersection.
fn intersect_segs(a0: V2, a1: V2, b0: V2, b1: V2) -> V2 {
    let d = ((a0.x - a1.x) * (b0.y - b1.y)) - ((a0.y - a1.y) * (b0.x - b1.x));

    if d.abs() < 0.000001 {
        return v2(NAN, NAN);
    }

    let t = (((a0.x - b0.x) * (b0.y - b1.y)) - ((a0.y - b0.y) * (b0.x - b1.x))) / d;
    let u = (((a0.x - b0.x) * (a0.y - a1.y)) - ((a0.y - b0.y) * (a0.x - a1.x))) / d;
    if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
        v2(a0.x + (t * (a1.x - a0.x)), a0.y + (t * (a1.y - a0.y)))
    } else {
        v2(NAN, NAN)
    }
}

fn abgr_mul(col: u32, a: u32) -> u32 {
    let br = ((col & 0xFF00FF) * a) >> 8;
    let g = ((col & 0x00FF00) * a) >> 8;
    0xFF000000 | (br & 0xFF00FF) | (g & 0x00FF00)
}

#[derive(Copy, Clone, Default, FromStr)]
#[display("{a_x} {a_y} {b_x} {b_y} {portal}")]
struct Wall {
    a_x: i32,
    a_y: i32,
    b_x: i32,
    b_y: i32,
    portal: usize,
}

// sector id for "no sector"
const SECTOR_NONE: usize = 0;
const SECTOR_MAX: usize = 128;

#[derive(Copy, Clone, Default, FromStr)]
#[display("{id} {firstwall} {nwalls} {zfloor} {zceil}")]
struct Sector {
    id: i32,
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

#[derive(Default)]
struct Camera {
    pos: V2,
    angle: f32,
    anglecos: f32,
    anglesin: f32,
    sector: usize,
}

struct State {
    pixels: Vec<u32>,
    level: Level,
    camera: Camera,
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

// world space -> camera space (translate and rotate)
fn world_pos_to_camera(state: &State, p: V2) -> V2 {
    let u = v2(p.x - state.camera.pos.x, p.y - state.camera.pos.y);
    v2(
        u.x * state.camera.anglesin - u.y * state.camera.anglecos,
        u.x * state.camera.anglecos + u.y * state.camera.anglesin,
    )
}

// load sectors from file -> state
fn load_level(path: &str) -> Result<Level, String> {
    let mut level = Level::default();

    // sector 0 does not exist
    level.sectors_count = 1;

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

fn verline(pixels: &mut [u32], x: usize, y0: i32, y1: i32, abgr: u32) {
    for y in (y0 as usize)..=(y1 as usize) {
        pixels[(SCREEN_HEIGHT - y - 1) * SCREEN_WIDTH + x] = abgr.to_be().rotate_right(8);
    }
}

// point is in sector if it is on the left side of all walls
fn point_in_sector(level: &Level, sector: &Sector, p: V2) -> bool {
    for i in 0..sector.nwalls {
        let wall = &level.walls[sector.firstwall + i];

        if point_side(
            p,
            v2(wall.a_x as f32, wall.a_y as f32),
            v2(wall.b_x as f32, wall.b_y as f32),
        ) > 0.0
        {
            return false;
        }
    }

    true
}

#[derive(Copy, Clone, Default)]
struct QueueEntry {
    id: usize,
    x0: usize,
    x1: usize,
}

fn render(state: &mut State) {
    let mut y_lo = [0u16; SCREEN_WIDTH];
    let mut y_hi = [(SCREEN_HEIGHT - 1) as u16; SCREEN_WIDTH];

    // track if sector has already been drawn
    let mut sectdraw = [false; SECTOR_MAX];

    // calculate edges of near/far planes (looking down +Y axis)
    let zdl = rotate(v2(0.0, 1.0), FRAC_HFOV_2);
    let zdr = rotate(v2(0.0, 1.0), -FRAC_HFOV_2);
    let znl = v2(zdl.x * ZNEAR, zdl.y * ZNEAR);
    let znr = v2(zdr.x * ZNEAR, zdr.y * ZNEAR);
    let zfl = v2(zdl.x * ZFAR, zdl.y * ZFAR);
    let zfr = v2(zdr.x * ZFAR, zdr.y * ZFAR);

    let mut queue = vec![QueueEntry {
        id: state.camera.sector,
        x0: 0,
        x1: SCREEN_WIDTH - 1,
    }];

    while !queue.is_empty() {
        // grab tail of queue
        let entry = queue.pop().unwrap();

        if sectdraw[entry.id] {
            continue;
        }

        sectdraw[entry.id] = true;

        let sector = &state.level.sectors[entry.id];

        for i in 0..sector.nwalls {
            let wall = &state.level.walls[sector.firstwall + i];

            // translate relative to player and rotate points around player's view
            let op0 = world_pos_to_camera(state, v2(wall.a_x as f32, wall.a_y as f32));
            let op1 = world_pos_to_camera(state, v2(wall.b_x as f32, wall.b_y as f32));

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
                if !il.x.is_nan() {
                    cp0 = il;
                    ap0 = normalize_angle(f32::atan2(cp0.y, cp0.x) - FRAC_PI_2);
                }

                if !ir.x.is_nan() {
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
            if tx0 > entry.x1 {
                continue;
            }
            if tx1 < entry.x0 {
                continue;
            }

            let wallshade =
                16.0 * f32::sin(f32::atan2(
                    (wall.b_x - wall.a_x) as f32,
                    (wall.b_y - wall.a_y) as f32,
                )) + 1.0;

            let x0 = usize::clamp(tx0, entry.x0, entry.x1);
            let x1 = usize::clamp(tx1, entry.x0, entry.x1);

            let z_floor = sector.zfloor;
            let z_ceil = sector.zceil;
            let nz_floor = if wall.portal > 0 {
                state.level.sectors[wall.portal].zfloor
            } else {
                0.0
            };
            let nz_ceil = if wall.portal > 0 {
                state.level.sectors[wall.portal].zceil
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
                let yf = i32::clamp(tyf, y_lo[x] as i32, y_hi[x] as i32);
                let yc = i32::clamp(tyc, y_lo[x] as i32, y_hi[x] as i32);

                // floor
                if yf > y_lo[x] as i32 {
                    verline(&mut state.pixels, x, y_lo[x] as i32, yf, 0xFFFF0000);
                }

                // ceiling
                if yc < y_hi[x] as i32 {
                    verline(&mut state.pixels, x, yc, y_hi[x] as i32, 0xFF00FFFF);
                }

                if wall.portal > 0 {
                    let tnyf = (xp * nyfd as f32) as i32 + nyf0;
                    let tnyc = (xp * nycd as f32) as i32 + nyc0;
                    let nyf = i32::clamp(tnyf, y_lo[x] as i32, y_hi[x] as i32);
                    let nyc = i32::clamp(tnyc, y_lo[x] as i32, y_hi[x] as i32);

                    verline(&mut state.pixels, x, nyc, yc, abgr_mul(0xFF00FF00, shade));
                    verline(&mut state.pixels, x, yf, nyf, abgr_mul(0xFF0000FF, shade));

                    y_hi[x] =
                        i32::clamp(yc.min(nyc).min(y_hi[x] as i32), 0, SCREEN_HEIGHT as i32 - 1)
                            as u16;
                    y_lo[x] =
                        i32::clamp(yf.max(nyf).max(y_lo[x] as i32), 0, SCREEN_HEIGHT as i32 - 1)
                            as u16;
                } else {
                    verline(&mut state.pixels, x, yf, yc, abgr_mul(0xFFD0D0D0, shade));
                }
            }

            if wall.portal > 0 {
                queue.push(QueueEntry {
                    id: wall.portal,
                    x0,
                    x1,
                });
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

    let pixels = vec![0; SCREEN_WIDTH * SCREEN_HEIGHT * 4];

    let level = load_level("level.txt").unwrap();

    let mut state = State {
        pixels,
        camera: Camera {
            pos: v2(3.0, 3.0),
            angle: 0.0,
            sector: 1,
            ..Default::default()
        },
        level,
    };

    while window.is_open() && !window.is_key_down(Key::Escape) {
        let rot_speed = 3.0 * 0.016;
        let move_speed = 3.0 * 0.016;

        if window.is_key_down(Key::Right) {
            state.camera.angle -= rot_speed;
        }

        if window.is_key_down(Key::Left) {
            state.camera.angle += rot_speed;
        }

        state.camera.anglecos = f32::cos(state.camera.angle);
        state.camera.anglesin = f32::sin(state.camera.angle);

        if window.is_key_down(Key::Up) {
            state.camera.pos = v2(
                state.camera.pos.x + (move_speed * state.camera.anglecos),
                state.camera.pos.y + (move_speed * state.camera.anglesin),
            );
        }

        if window.is_key_down(Key::Down) {
            state.camera.pos = v2(
                state.camera.pos.x - (move_speed * state.camera.anglecos),
                state.camera.pos.y - (move_speed * state.camera.anglesin),
            );
        }

        // update player sector
        {
            // BFS neighbors in a circular queue, player is likely to be in one
            // of the neighboring sectors

            let mut queue = VecDeque::new();
            queue.push_back(state.camera.sector);
            let mut found = SECTOR_NONE;

            while let Some(id) = queue.pop_front() {
                let sector = &state.level.sectors[id];

                if point_in_sector(&state.level, sector, state.camera.pos) {
                    found = id;
                    break;
                }

                // check neighbors
                for j in 0..sector.nwalls {
                    let wall = &state.level.walls[sector.firstwall + j];

                    if wall.portal > 0 {
                        queue.push_back(wall.portal)
                    }
                }
            }

            if found == SECTOR_NONE {
                println!("Player is not in a sector!");
                state.camera.sector = 1;
            } else {
                state.camera.sector = found;
            }
        }

        render(&mut state);
        window
            .update_with_buffer(&state.pixels, SCREEN_WIDTH, SCREEN_HEIGHT)
            .unwrap();
    }
}
