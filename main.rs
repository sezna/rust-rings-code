extern crate rand;
struct vec3d {
    x:f32,
    y:f32,
    z:f32,
}
struct ray {
    r0: vec3d,
    r: vec3d,
}

struct particle {
    x: f32,
    y: f32,
    z: f32,
    vx: f32,
    vy: f32,
    vz: f32,
    rad: f32,
}

struct bindata {
    xmin: f32,
    xmax: f32,
    ymin: f32,
    ymax: f32,
    bin_size: f32,
    bins: Vec<Vec<Vec<particle>>>,
}

struct photon {
    x: f32,
    y: f32,
    hit: bool,
}

struct scan {
    sx: f32,
    sy: f32,
    ex: f32,
    ey: f32,
    intensity: f32,
    photons: Vec<photon>,
}


fn main() {
    
}

fn synthetic_occultation(x: f32, y: f32, theta: f32, phi: f32, cut_theta: f32, scan_length: f32, off_length: f32, beam_size: f32, bin_data: bindata, zmin: f32, zmax: f32, photon_count: i32) -> Vec<scan> {
    let r_dir:vec3d = vec3d {
        x: theta.cos() * phi.cos(), 
        y: theta.sin() * phi.cos(),
        z: phi.sin(),
    };
    let dx = cut_theta.cos();
    let dy = cut_theta.sin();
    let mut height = 0 as f32;
    if zmin.abs() > zmax {
     height = zmin.abs();   
    }
    else {
     height = zmax;   
    }
    let xstart = bin_data.xmin + rDir.x * height;
    let xend = bin_data.xmax - rDir * height;
    let mx = xstart;
    let ret:Vec<scan> = Vec::new();
    let mut sx:f32;
    let mut sy:f32;
    let mut ex:f32;
    let mut ey:f32;
    let mut photons:Vec<photon> = Vec::new();
    let mut t:f32;
    let mut rx:f32;
    let mut ry:f32;
    fn photon_map_func(x: i32) -> photon{
        let t = rand::random::<f64>();
        let rx = sx + t * (ex - sx) + rand::random::<f64>() * rand::random::<f64>() * beam_size;
        let ry = sy + t * (ey - sy) + rand::random::<f64>() * rand::random::<f64>() * beam_size;
        return photon {x: rx, 
                       y: ry, 
                       hit: ray_grid_intersect(ray{r0: vec3d{x: rx, y: ry, z: 0}, r: r_dir}, bin_data, zmin, zmax),
        };
    }
    while mx < xend {
        sx = mx;
        sy = y + (cut_theta.tan() * (sx - x));
        ex = sx + scan_length * cut_theta.cos();
        ey = y + cut_theta.tan() * (ex - x);
        photons = (1..photon_count).map(photon_map_func(_)); 
    }
    ret.insert(scan{sx: sx, sy: sy, ex: ex, ey: ey, intensity: photons.count(p => !p.hit) as f32 / photon_count, photons:
    photons.seq});
    mx = mx + (scan_length + off_length) * cut_theta.cos();
    return ret;
}

fn ray_grid_intersect(r: ray, bin_data: bindata, zmin: f32, zmax: f32) -> bool {

    return false;
}

fn ray_bin_intersect(r: ray, xbin: i32, ybin: i32, bin_data: bindata) -> bool {

    return false;
}

fn ray_particle_intersect(r: ray, part: particle) -> bool {

    return false;
}

fn bin_particles (parts: Vec<particle>) -> bindata {
    
    let tmpvec:Vec<Vec<Vec<particle>>> = Vec::new();
    return bindata {
    xmin: 0 as f32,
    xmax: 0 as f32,
    ymin: 0 as f32,
    ymax: 0 as f32,
    bin_size: 0 as f32,
    bins: tmpvec,
    }
}


