use std::io::Cursor;
use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::Read;
extern crate byteorder;
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
#[derive(Clone)]
struct cart {
    x: f32,
    y: f32,
    z: f32,
    vx: f32,
    vy: f32,
    vz: f32,
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

impl ray {
    fn clone(&self) -> ray {
        return ray {
            r0: self.r0.clone(),
            r: self.r.clone(),
        }
    }
}
impl vec3d {
    fn times(&self, other: f32) -> vec3d {
        return vec3d {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
    fn clone(&self) -> vec3d {
        return vec3d {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
    fn sub(&self, other: &vec3d) -> vec3d {
        return vec3d {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    fn dot(&self, other: &vec3d) -> f32 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }
}

fn main() {
  let xbindensity = 10000;
  let ybindensity = 100;
  let file_name = "test_output.bin";
    //read in file
    let mut file = match File::open(file_name) {
        Ok(file) => file,
        Err(..) => {
            println!("could not open file {}", file_name);
            return;
        },
    };
    let num:i32 = file.read_i32::<LittleEndian>().unwrap();
    println!("There are {} particles", num);
    let mut tmpcart = cart {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        vx: 0.0,
        vy: 0.0,
        vz: 0.0,
    };
    let mut cartesians:Vec<cart> = Vec::new();
    for i in (1..num) {
        let mut tmpcart = cart {
               x: 0.0,
                y: 0.0,
                z: 0.0,
                vx: 0.0,
                vy: 0.0,
                vz: 0.0,
            };
         tmpcart.x = file.read_f32::<LittleEndian>().unwrap();
         tmpcart.y = file.read_f32::<LittleEndian>().unwrap();
         tmpcart.z = file.read_f32::<LittleEndian>().unwrap();
         tmpcart.vx = file.read_f32::<LittleEndian>().unwrap();
         tmpcart.vy = file.read_f32::<LittleEndian>().unwrap();
         tmpcart.vz = file.read_f32::<LittleEndian>().unwrap();
         cartesians.push(tmpcart);
    }
    let mut rads:Vec<f32> = Vec::new();
    let mut tmprad:f32 = 0.0;
    for i in (1..num) {
        tmprad = file.read_f32::<LittleEndian>().unwrap();
        rads.push(tmprad);
    }
//sizeof int: 4 sizeof cart: 48, sizeof double: 8x
    println!("first cart: x: {} y: {} z: {} vx: {} vy: {} vz: {}", cartesians[0].x,
             cartesians[0].y, cartesians[0].z, cartesians[0].vx, cartesians[0].vy,
             cartesians[0].vz);
// Find xmin, xmax, ymin, ymax
    let mut xmax = cartesians[0].x;
    let mut xmin = cartesians[0].x;
    let mut ymax = cartesians[0].y;
    let mut ymin = cartesians[0].y;
  for i in (0..num as usize - 1 as usize) {
    if cartesians[i].x < xmin {
        xmin = cartesians[i].x;
    }
    if cartesians[i].y < ymin {
        ymin = cartesians[i].y;
    }
    if cartesians[i].x > xmax {
        xmax = cartesians[i].x;
    }
    if cartesians[i].y > ymax {
        ymax = cartesians[i].y
    }
  }
  println!("xmax: {} ymax: {} xmin: {} ymin: {}", xmax, ymax, xmin, ymin);
    let mut grid:Vec<Vec<Vec<i32>>> = vec![vec![vec![]; ybindensity]; xbindensity];
   for i in (0..num as usize - 1 as usize) {
        let xlen = grid.len();
        let ylen = grid[0].len();
        let mut xbin:i32 = ((cartesians[i].x - xmin) * grid.len() as f32 / (xmax - xmin)) as i32;
        let mut ybin:i32 = ((cartesians[i].y - ymin) * grid[0].len() as f32 / (ymax - ymin)) as i32;
        if xbin as usize >= xlen {
            xbin  = xlen as i32 - 1;
        }
        if ybin as usize >= ylen {
            ybin  = ylen as i32 - 1;
        }
        let xylen = grid[xbin as usize][ybin as usize].len();
        grid[xbin as usize][ybin as usize].insert(xylen, i as i32);
   }

}


fn synthetic_occultation(x: f32, y: f32, theta: f32, phi: f32, cut_theta: f32, scan_length: f32, off_length: f32, beam_size: f32, bin_data: &bindata, zmin: f32, zmax: f32, photon_count: i32) -> Vec<scan> {
    let  r_dir:vec3d = vec3d {
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
    let xstart = bin_data.xmin + r_dir.x * height;
    let xend = bin_data.xmax - r_dir.x * height;
    let mut mx = xstart;
    let mut ret:Vec<scan> = Vec::new();
    let mut sx:f32 = 0.0;
    let mut sy:f32 = 0.0;
    let mut ex:f32 = 0.0;
    let mut ey:f32 = 0.0;
    let mut photons:Vec<photon> = Vec::new();
    let mut t:f32 = 0.0;
    let mut rx:f32 = 0.0;
    let mut ry:f32 = 0.0;
    while mx < xend {
        sx = mx;
        sy = y + (cut_theta.tan() * (sx - x));
        ex = sx + scan_length * cut_theta.cos();
        ey = y + cut_theta.tan() * (ex - x);
        for i in (1..photon_count) {
            let t = rand::random::<f32>();
            let rx = sx + t * (ex - sx) + rand::random::<f32>() * rand::random::<f32>() * beam_size;
            let ry = sy + t * (ey - sy) + rand::random::<f32>() * rand::random::<f32>() * beam_size;
            photons.push(photon {x: rx, y: ry, hit: ray_grid_intersect(ray{r0: vec3d{x: rx, y: ry, z: 0 as f32}, r: r_dir.clone()}, bin_data, zmin, zmax)})

        }
    }
    ret.push(scan{sx: sx, sy: sy, ex: ex, ey: ey, intensity: photons.iter_mut().filter(|p|{!p.hit}).count() as f32 / photon_count as f32, photons:
    photons});
    mx = mx + (scan_length + off_length) * cut_theta.cos();
    return ret;
}

fn ray_grid_intersect(r: ray, bin_data: &bindata, zmin: f32, zmax: f32) -> bool {
    let tmin = (zmin - r.r0.z) / r.r.z;
    let tmax = (zmax - r.r0.z) / r.r.z;
    let xmin = r.r0.x + tmin * r.r.x;
    let xmax = r.r0.x + tmax * r.r.x;
    let ymin = r.r0.y + tmin * r.r.y;
    let ymax = r.r0.y + tmax * r.r.y;
    let minxbin = max(((min(xmin, xmax) - bin_data.xmin) / bin_data.bin_size) - 1 as f32, 0 as f32)
        as i32;
    let maxxbin = min(((max(xmin, xmax) - bin_data.xmin) / bin_data.bin_size) + 1 as f32,
    bin_data.bins.len() as f32 - 1 as f32) as i32;
    let minybin = max(((min(ymin, ymax) - bin_data.ymin) / bin_data.bin_size) - 1 as f32, 0 as f32)
        as i32;
    let maxybin = min(((max(ymin, ymax) - bin_data.ymin) / bin_data.bin_size) + 1 as f32,
    bin_data.bins[0].len() as f32 - 1 as f32) as i32;
    for xbin in (minxbin..maxxbin) {
        for ybin in (minybin..maxybin) {
            if ray_bin_intersect(r.clone(), xbin, ybin, bin_data) {
                return true;   
            }
        }
    }
    return false;
}

fn ray_bin_intersect(r: ray, xbin: i32, ybin: i32, bin_data: &bindata) -> bool {
        for i in bin_data.bins[xbin as usize][ybin as usize].iter() {
            if (ray_particle_intersect(&r, i)) {
                return true;
            }
        }
    return false;
}

fn ray_particle_intersect(r: &ray, part: &particle) -> bool {
    let p = vec3d{
        x: part.x,
        y: part.y,
        z: part.z,
    };
    let d = vec3d::sub(&r.r0, &p);
    let a = vec3d::dot(&r.r, &r.r);
    let b = vec3d::dot(&r.r, &d) * 2 as f32;
    let c = vec3d::dot(&d, &d) - part.rad * part.rad;
    return b * b - 4 as f32 * a * c >= 0 as f32;
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

fn min(a: f32, b:f32) -> f32 {
    if a < b {
        return a;
    }
    else {
        return b;
    }
}
fn max(a: f32, b:f32) -> f32 {
    if a > b {
        return a;
    }
    else {
        return b;
    }
}
