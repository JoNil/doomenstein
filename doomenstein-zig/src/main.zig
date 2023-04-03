const std = @import("std");
const math = @import("std").math;

const SCREEN_WIDTH: usize = 384;
const SCREEN_HEIGHT: usize = 216;

const EYE_Z: f32 = 1.65;
const HFOV: f32 = 120.0 * math.pi / 180.0;
const FRAC_HFOV_2: f32 = HFOV / 2.0;
const VFOV: f32 = 0.5;

const ZNEAR: f32 = 0.0001;
const ZFAR: f32 = 128.0;

const SECTOR_NONE: usize = 0;
const SECTOR_MAX: usize = 128;

const Vec2 = struct {
    x: f32,
    y: f32,
};

const Mat2 = struct {
    m00: f32,
    m01: f32,
    m10: f32,
    m11: f32,
};

const Camera = struct {
    pos: Vec2,
    angle: f32,
    camera_rot: Mat2,
    sector: usize,
};

const Wall = struct {
    a: Vec2,
    b: Vec2,
    portal: usize,
};

const Sector = struct {
    _id: i32,
    firstwall: usize,
    nwalls: usize,
    zfloor: f32,
    zceil: f32,
};

const Level = struct {
    sectors: [32]Sector,
    sectors_count: usize,
    walls: [32]Wall,
    walls_count: usize,
};

fn ifnan(x: f32, alt: f32) f32 {
    if (math.isNan(x)) return alt;
    return x;
}

fn intersect_segs(a1: Vec2, a2: Vec2, b1: Vec2, b2: Vec2) ?Vec2 {
    const u = a2 - a1;
    const v = b2 - b1;
    const w = a1 - b1;

    const d = (u.x * v.y) - (u.y * v.x);

    if (d == 0.0) return null;

    const t = ((v.y * w.x) - (v.x * w.y)) / d;
    const s = ((u.y * w.x) - (u.x * w.y)) / d;

    if (t < 0.0 or t > 1.0 or s < 0.0 or s > 1.0) return null;

    return Vec2{
        .x = a1.x + (t * u.x),
        .y = a1.y + (t * u.y),
    };
}

fn abgr_mul(a: u32, b: u32) u32 {
    const a_a = (a >> 24) & 0xFF;
    const a_r = (a >> 16) & 0xFF;
    const a_g = (a >> 8) & 0xFF;
    const a_b = a & 0xFF;

    const b_a = (b >> 24) & 0xFF;
    const b_r = (b >> 16) & 0xFF;
    const b_g = (b >> 8) & 0xFF;
    const b_b = b & 0xFF;

    const r_a = a_a * b_a / 255;
    const r_r = a_r * b_r / 255;
    const r_g = a_g * b_g / 255;
    const r_b = a_b * b_b / 255;

    return (r_a << 24) | (r_r << 16) | (r_g << 8) | r_b;
}

pub fn main() !void {
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});
}
