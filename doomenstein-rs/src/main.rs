use std::{
    f32::{
        consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU},
        NAN,
    },
    fs::{self, File},
    path::Path,
};

use parse_display::FromStr;
use sdl2::{
    render::{RendererContext, Texture},
    video::{Window, WindowContext},
};

const SCREEN_WIDTH: usize = 384;
const SCREEN_HEIGHT: usize = 216;

const EYE_Z: f32 = 1.65;
const HFOV: f32 = 90.0 * std::f32::consts::PI / 180.0;
const VFOV: f32 = 0.5;

const ZNEAR: f32 = 0.0001;
const ZFAR: f32 = 128.0;

#[derive(Copy, Clone, FromStr)]
#[display("{x} {y}")]
struct v2 {
    x: f32,
    y: f32,
}

#[derive(Copy, Clone, FromStr)]
#[display("{x} {y}")]
struct v2i {
    x: i32,
    y: i32,
}

fn v2_to_v2i(v: v2) -> v2i {
    v2i {
        x: v.x as i32,
        y: v.y as i32,
    }
}

fn v2i_to_v2(v: v2i) -> v2 {
    v2 {
        x: v.x as f32,
        y: v.y as f32,
    }
}

fn dot(v0: v2, v1: v2) -> f32 {
    (v0.x * v1.x) + (v0.y * v1.y)
}

fn length(v: v2) -> f32 {
    dot(v, v).sqrt()
}

fn noramlize(v: v2) -> v2 {
    let len = length(v);
    v2 {
        x: v.x / len,
        y: v.y / len,
    }
}

fn ifnan(x: f32, alt: f32) -> f32 {
    if x.is_nan() {
        alt
    } else {
        x
    }
}

// -1 right, 0 on, 1 left
fn point_side(p: v2, a: v2, b: v2) -> f32 {
    -(((p.x - a.x) * (b.y - a.y)) - ((p.y - a.y) * (b.x - a.x)))
}

// rotate vector v by angle a
fn rotate(v: v2, a: f32) -> v2 {
    v2 {
        x: (v.x * a.cos()) - (v.y * a.sin()),
        y: (v.x * a.sin()) + (v.y * a.cos()),
    }
}

// see: https://en.wikipedia.org/wiki/Line–line_intersection
// compute intersection of two line segments, returns (NAN, NAN) if there is
// no intersection.
fn intersect_segs(a0: v2, a1: v2, b0: v2, b1: v2) -> v2 {
    let d = ((a0.x - a1.x) * (b0.y - b1.y)) - ((a0.y - a1.y) * (b0.x - b1.x));

    if d.abs() < 0.000001 {
        return v2 { x: NAN, y: NAN };
    }

    let t = (((a0.x - b0.x) * (b0.y - b1.y)) - ((a0.y - b0.y) * (b0.x - b1.x))) / d;
    let u = (((a0.x - b0.x) * (a0.y - a1.y)) - ((a0.y - b0.y) * (a0.x - a1.x))) / d;
    if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
        v2 {
            x: a0.x + (t * (a1.x - a0.x)),
            y: a0.y + (t * (a1.y - a0.y)),
        }
    } else {
        v2 { x: NAN, y: NAN }
    }
}

fn abgr_mul(col: u32, a: u32) -> u32 {
    let br = ((col & 0xFF00FF) * a) >> 8;
    let g = ((col & 0x00FF00) * a) >> 8;

    0xFF000000 | (br & 0xFF00FF) | (g & 0x00FF00)
}

#[derive(FromStr)]
#[display("{a} {b} {portal}")]
struct Wall {
    a: v2i,
    b: v2i,
    portal: i32,
}

// sector id for "no sector"
const SECTOR_NONE: u32 = 0;
const SECTOR_MAX: u32 = 128;

#[derive(FromStr)]
#[display("{id} {firstwall} {nwalls} {zfloor} {zceil}")]
struct Sector {
    id: i32,
    firstwall: usize,
    nwalls: usize,
    zfloor: f32,
    zceil: f32,
}

struct Camera {
    pos: v2,
    angle: f32,
    anglecos: f32,
    anglesin: f32,
    sector: i32,
}

struct State {
    window: Window,
    renderer: RendererContext<WindowContext>,
    texture: Texture,
    debug: Texture,
    pixels: Vec<u32>,
    quit: bool,

    sectors: [Sector; 32],
    sectors_count: usize,

    walls: [Wall; 32],
    walls_count: usize,

    y_lo: [u16; SCREEN_WIDTH],
    y_hi: [u16; SCREEN_WIDTH],

    camera: Camera,

    sleepy: bool,
}

// convert angle in [-(HFOV / 2)..+(HFOV / 2)] to X coordinate
fn screen_angle_to_x(angle: f32) -> i32 {
    ((SCREEN_WIDTH as f32 / 2.0)
        * (1.0 - (((angle + (HFOV / 2.0)) / HFOV) * FRAC_PI_2 - FRAC_PI_4).tan())) as i32
}

// noramlize angle to +/-PI
fn normalize_angle(a: f32) -> f32 {
    return a - (TAU * ((a + PI) / TAU).floor());
}

// world space -> camera space (translate and rotate)
fn world_pos_to_camera(state: &State, p: v2) -> v2 {
    let u = v2 {
        x: p.x - state.camera.pos.x,
        y: p.y - state.camera.pos.y,
    };
    v2 {
        x: u.x * state.camera.anglesin - u.y * state.camera.anglecos,
        y: u.x * state.camera.anglecos + u.y * state.camera.anglesin,
    }
}

// load sectors from file -> state
fn load_sectors(state: &mut State, path: &Path) -> Result<(), String> {
    // sector 0 does not exist
    state.sectors_count = 1;

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
        }
        if line == "[WALL]" {
            parse_state = State::Wall;
        }

        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        match parse_state {
            State::Sector => {
                let sector = line
                    .parse()
                    .map_err(|e| format!("Unable to parse sector {path:?}:{i}: {e}"))?;
                state.sectors[state.sectors_count] = sector;
                state.sectors_count += 1;
            }
            State::Wall => {
                let wall = line
                    .parse()
                    .map_err(|e| format!("Unable to parse wall {path:?}:{i}: {e}"))?;
                state.walls[state.walls_count] = wall;
                state.walls_count += 1;
            }
        }
    }

    Ok(())
}

fn verline(state: &mut State, x: i32, y0: i32, y1: i32, color: u32) {
    for y in y0..=y1 {
        state.pixels[(y * SCREEN_WIDTH as i32 + x) as usize] = color;
    }
}

/*
// point is in sector if it is on the left side of all walls
static bool point_in_sector(const struct sector *sector, v2 p) {
    for (usize i = 0; i < sector->nwalls; i++) {
        const struct wall *Wall = &state.walls.arr[sector->firstwall + i];

        if (point_side(p, v2i_to_v2(Wall->a), v2i_to_v2(Wall->b)) > 0) {
            return false;
        }
    }

    return true;
}

static void render() {
    for (int i = 0; i < SCREEN_WIDTH; i++) {
        state.y_hi[i] = SCREEN_HEIGHT - 1;
        state.y_lo[i] = 0;
    }

    // track if sector has already been drawn
    bool sectdraw[SECTOR_MAX];
    memset(sectdraw, 0, sizeof(sectdraw));

    // calculate edges of near/far planes (looking down +Y axis)
    const v2
        zdl = rotate(((v2) { 0.0f, 1.0f }), +(HFOV / 2.0f)),
        zdr = rotate(((v2) { 0.0f, 1.0f }), -(HFOV / 2.0f)),
        znl = (v2) { zdl.x * ZNEAR, zdl.y * ZNEAR },
        znr = (v2) { zdr.x * ZNEAR, zdr.y * ZNEAR },
        zfl = (v2) { zdl.x * ZFAR, zdl.y * ZFAR },
        zfr = (v2) { zdr.x * ZFAR, zdr.y * ZFAR };

    enum { QUEUE_MAX = 64 };
    struct queue_entry { int id, x0, x1; };

    struct { struct queue_entry arr[QUEUE_MAX]; usize n; } queue = {
        {{ state.camera.sector, 0, SCREEN_WIDTH - 1 }},
        1
    };

    while (queue.n != 0) {
        // grab tail of queue
        struct queue_entry entry = queue.arr[--queue.n];

        if (sectdraw[entry.id]) {
            continue;
        }

        sectdraw[entry.id] = true;

        const struct sector *sector = &state.sectors.arr[entry.id];

        for (usize i = 0; i < sector->nwalls; i++) {
            const struct wall *Wall =
                &state.walls.arr[sector->firstwall + i];

            // translate relative to player and rotate points around player's view
            const v2
                op0 = world_pos_to_camera(v2i_to_v2(Wall->a)),
                op1 = world_pos_to_camera(v2i_to_v2(Wall->b));

            // wall clipped pos
            v2 cp0 = op0, cp1 = op1;

            // both are negative -> wall is entirely behind player
            if (cp0.y <= 0 && cp1.y <= 0) {
                continue;
            }

            // angle-clip against view frustum
            f32
                ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2),
                ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);

            // clip against view frustum if both angles are not clearly within
            // HFOV
            if (cp0.y < ZNEAR
                || cp1.y < ZNEAR
                || ap0 > +(HFOV / 2)
                || ap1 < -(HFOV / 2)) {
                const v2
                    il = intersect_segs(cp0, cp1, znl, zfl),
                    ir = intersect_segs(cp0, cp1, znr, zfr);

                // recompute angles if points change
                if (!isnan(il.x)) {
                    cp0 = il;
                    ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2);
                }

                if (!isnan(ir.x)) {
                    cp1 = ir;
                    ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);
                }
            }

            if (ap0 < ap1) {
                continue;
            }

            if ((ap0 < -(HFOV / 2) && ap1 < -(HFOV / 2))
                || (ap0 > +(HFOV / 2) && ap1 > +(HFOV / 2))) {
                continue;
            }

            // "true" xs before portal clamping
            const int
                tx0 = screen_angle_to_x(ap0),
                tx1 = screen_angle_to_x(ap1);

            // bounds check against portal window
            if (tx0 > entry.x1) { continue; }
            if (tx1 < entry.x0) { continue; }

            const int wallshade =
                16 * (sin(atan2f(
                    Wall->b.x - Wall->a.x,
                    Wall->b.y - Wall->b.y)) + 1.0f);

            const int
                x0 = clamp(tx0, entry.x0, entry.x1),
                x1 = clamp(tx1, entry.x0, entry.x1);

            const f32
                z_floor = sector->zfloor,
                z_ceil = sector->zceil,
                nz_floor =
                    Wall->portal ? state.sectors.arr[Wall->portal].zfloor : 0,
                nz_ceil =
                    Wall->portal ? state.sectors.arr[Wall->portal].zceil : 0;

            const f32
                sy0 = ifnan((VFOV * SCREEN_HEIGHT) / cp0.y, 1e10),
                sy1 = ifnan((VFOV * SCREEN_HEIGHT) / cp1.y, 1e10);

            const int
                yf0  = (SCREEN_HEIGHT / 2) + (int) (( z_floor - EYE_Z) * sy0),
                yc0  = (SCREEN_HEIGHT / 2) + (int) (( z_ceil  - EYE_Z) * sy0),
                yf1  = (SCREEN_HEIGHT / 2) + (int) (( z_floor - EYE_Z) * sy1),
                yc1  = (SCREEN_HEIGHT / 2) + (int) (( z_ceil  - EYE_Z) * sy1),
                nyf0 = (SCREEN_HEIGHT / 2) + (int) ((nz_floor - EYE_Z) * sy0),
                nyc0 = (SCREEN_HEIGHT / 2) + (int) ((nz_ceil  - EYE_Z) * sy0),
                nyf1 = (SCREEN_HEIGHT / 2) + (int) ((nz_floor - EYE_Z) * sy1),
                nyc1 = (SCREEN_HEIGHT / 2) + (int) ((nz_ceil  - EYE_Z) * sy1),
                txd = tx1 - tx0,
                yfd = yf1 - yf0,
                ycd = yc1 - yc0,
                nyfd = nyf1 - nyf0,
                nycd = nyc1 - nyc0;

            for (int x = x0; x <= x1; x++) {
                int shade = x == x0 || x == x1 ? 192 : (255 - wallshade);

                // calculate progress along x-axis via tx{0,1} so that walls
                // which are partially cut off due to portal edges still have
                // proper heights
                const f32 xp = ifnan((x - tx0) / (f32) txd, 0);

                // get y coordinates for this x
                const int
                    tyf = (int) (xp * yfd) + yf0,
                    tyc = (int) (xp * ycd) + yc0,
                    yf = clamp(tyf, state.y_lo[x], state.y_hi[x]),
                    yc = clamp(tyc, state.y_lo[x], state.y_hi[x]);

                // floor
                if (yf > state.y_lo[x]) {
                    verline(
                        x,
                        state.y_lo[x],
                        yf,
                        0xFFFF0000);
                }

                // ceiling
                if (yc < state.y_hi[x]) {
                    verline(
                        x,
                        yc,
                        state.y_hi[x],
                        0xFF00FFFF);
                }

                if (Wall->portal) {
                    const int
                        tnyf = (int) (xp * nyfd) + nyf0,
                        tnyc = (int) (xp * nycd) + nyc0,
                        nyf = clamp(tnyf, state.y_lo[x], state.y_hi[x]),
                        nyc = clamp(tnyc, state.y_lo[x], state.y_hi[x]);

                    verline(
                        x,
                        nyc,
                        yc,
                        abgr_mul(0xFF00FF00, shade));

                    verline(
                        x,
                        yf,
                        nyf,
                        abgr_mul(0xFF0000FF, shade));

                    state.y_hi[x] =
                        clamp(
                            min(min(yc, nyc), state.y_hi[x]),
                            0, SCREEN_HEIGHT - 1);

                    state.y_lo[x] =
                        clamp(
                            max(max(yf, nyf), state.y_lo[x]),
                            0, SCREEN_HEIGHT - 1);
                } else {
                    verline(
                        x,
                        yf,
                        yc,
                        abgr_mul(0xFFD0D0D0, shade));
                }

                if (state.sleepy) {
                    present();
                    SDL_Delay(10);
                }
            }

            if (Wall->portal) {
                ASSERT(queue.n != QUEUE_MAX, "out of queue space");
                queue.arr[queue.n++] = (struct queue_entry) {
                    .id = Wall->portal,
                    .x0 = x0,
                    .x1 = x1
                };
            }
        }
    }

    state.sleepy = false;
}

static void present() {
    void *px;
    int pitch;
    SDL_LockTexture(state.texture, NULL, &px, &pitch);
    {
        for (usize y = 0; y < SCREEN_HEIGHT; y++) {
            memcpy(
                &((u8*) px)[y * pitch],
                &state.pixels[y * SCREEN_WIDTH],
                SCREEN_WIDTH * 4);
        }
    }
    SDL_UnlockTexture(state.texture);

    SDL_SetRenderTarget(state.renderer, NULL);
    SDL_SetRenderDrawColor(state.renderer, 0, 0, 0, 0xFF);
    SDL_SetRenderDrawBlendMode(state.renderer, SDL_BLENDMODE_NONE);

    SDL_RenderClear(state.renderer);
    SDL_RenderCopyEx(
        state.renderer,
        state.texture,
        NULL,
        NULL,
        0.0,
        NULL,
        SDL_FLIP_VERTICAL);

    SDL_SetTextureBlendMode(state.debug, SDL_BLENDMODE_BLEND);
    SDL_RenderCopy(state.renderer, state.debug, NULL, &((SDL_Rect) { 0, 0, 512, 512 }));
    SDL_RenderPresent(state.renderer);
}

int main(int argc, char *argv[]) {
    ASSERT(
        !SDL_Init(SDL_INIT_VIDEO),
        "SDL failed to initialize: %s",
        SDL_GetError());

    state.window =
        SDL_CreateWindow(
            "raycast",
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            1280,
            720,
            0);

    ASSERT(state.window, "failed to create SDL window: %s\n", SDL_GetError());

    state.renderer =
        SDL_CreateRenderer(
            state.window,
            -1,
            SDL_RENDERER_ACCELERATED
            | SDL_RENDERER_PRESENTVSYNC);

    state.texture =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_STREAMING,
            SCREEN_WIDTH,
            SCREEN_HEIGHT);
    state.debug =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_TARGET,
            128,
            128);

    state.pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);

    state.camera.pos = (v2) { 3, 3 };
    state.camera.angle = 0.0;
    state.camera.sector = 1;

    int ret = 0;
    ASSERT(
        !(ret = load_sectors("res/level.txt")),
        "error while loading sectors: %d",
        ret);
    printf(
        "loaded %zu sectors with %zu walls",
        state.sectors.n,
        state.walls.n);

    while (!state.quit) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
                case SDL_QUIT:
                    state.quit = true;
                    break;
                default:
                    break;
            }
        }

        if (state.quit) {
            break;
        }

        const f32 rot_speed = 3.0f * 0.016f, move_speed = 3.0f * 0.016f;

        const u8 *keystate = SDL_GetKeyboardState(NULL);

        if (keystate[SDLK_RIGHT & 0xFFFF]) {
            state.camera.angle -= rot_speed;
        }

        if (keystate[SDLK_LEFT & 0xFFFF]) {
            state.camera.angle += rot_speed;
        }

        state.camera.anglecos = cos(state.camera.angle);
        state.camera.anglesin = sin(state.camera.angle);

        if (keystate[SDLK_UP & 0xFFFF]) {
            state.camera.pos = (v2) {
                state.camera.pos.x + (move_speed * state.camera.anglecos),
                state.camera.pos.y + (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_DOWN & 0xFFFF]) {
            state.camera.pos = (v2) {
                state.camera.pos.x - (move_speed * state.camera.anglecos),
                state.camera.pos.y - (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_F1 & 0xFFFF]) {
            state.sleepy = true;
        }

        // update player sector
        {
            // BFS neighbors in a circular queue, player is likely to be in one
            // of the neighboring sectors
            enum { QUEUE_MAX = 64 };
            int
                queue[QUEUE_MAX] = { state.camera.sector },
                i = 0,
                n = 1,
                found = SECTOR_NONE;

            while (n != 0) {
                // get front of queue and advance to next
                const int id = queue[i];
                i = (i + 1) % (QUEUE_MAX);
                n--;

                const struct sector *sector = &state.sectors.arr[id];

                if (point_in_sector(sector, state.camera.pos)) {
                    found = id;
                    break;
                }

                // check neighbors
                for (usize j = 0; j < sector->nwalls; j++) {
                    const struct wall *Wall =
                        &state.walls.arr[sector->firstwall + j];

                    if (Wall->portal) {
                        if (n == QUEUE_MAX) {
                            fprintf(stderr, "out of queue space!");
                            goto done;
                        }

                        queue[(i + n) % QUEUE_MAX] = Wall->portal;
                        n++;
                    }
                }
            }


done:
            if (!found) {
                fprintf(stderr, "player is not in a sector!");
                state.camera.sector = 1;
            } else {
                state.camera.sector = found;
            }
        }

        memset(state.pixels, 0, SCREEN_WIDTH * SCREEN_HEIGHT * 4);
        render();
        if (!state.sleepy) { present(); }
    }

    SDL_DestroyTexture(state.debug);
    SDL_DestroyTexture(state.texture);
    SDL_DestroyRenderer(state.renderer);
    SDL_DestroyWindow(state.window);
    return 0;
}
*/
