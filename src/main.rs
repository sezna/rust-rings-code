/*TODO: handle command line arguments, run the actual occultations*/
use std::cmp::Ordering;
use std::io::Cursor;
use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::Read;
extern crate byteorder;
extern crate rand;
struct vec3D {
    x:f64,
    y:f64,
    z:f64,
}
struct ray {
    r0: vec3D,
    r: vec3D,
}
#[derive(Clone)]
struct cart {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

struct particle {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    rad: f64,
}
impl cart {
    fn to_particle(&self, radius: f64) -> particle{
        return particle {
            x: self.x,
            y: self.y,
            z: self.z,
            vx: self.vx,
            vz: self.vz,
            vy: self.vy,
            rad: radius,
        }
    }
}
struct binData {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    bin_size: f64,
    bins: Vec<Vec<Vec<i32>>>,
}

struct photon {
    x: f64,
    y: f64,
    hit: bool,
}

struct scan {
    sx: f64,
    sy: f64,
    ex: f64,
    ey: f64,
    intensity: f64,
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
impl vec3D {
    fn times(&self, other: f64) -> vec3D {
        return vec3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
    fn clone(&self) -> vec3D {
        return vec3D {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
    fn sub(&self, other: &vec3D) -> vec3D {
        return vec3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    fn dot(&self, other: &vec3D) -> f64 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }
}

fn main() {
//fn synthetic_occultation(x: f64, y: f64, theta: f64, phi: f64, cut_theta: f64, scan_length: f64, off_length: f64, beam_size: f64, bin_data: &binData, zmin: f64, zmax: f64, photon_count: i32) -> Vec<scan> {
  let x = 0.0;
  let y = 0.0;
  let theta = 0.0;
  let phi = 1.57;
  let cut_theta = 0.0;
  let scan_length = 100.0 / 136505500.0;
  let off_length = 20.0 / 136505500.0;
  let beam_size = 10.0 / 136505500.0;
  let photon_count = 1000;
  let xbindensity = 10000;
  let ybindensity = 100;
  let file_name = "CartAndRad.10000.bin";
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
         tmpcart.x = file.read_f64::<LittleEndian>().unwrap();
         tmpcart.y = file.read_f64::<LittleEndian>().unwrap();
         tmpcart.z = file.read_f64::<LittleEndian>().unwrap();
         tmpcart.vx = file.read_f64::<LittleEndian>().unwrap();
         tmpcart.vy = file.read_f64::<LittleEndian>().unwrap();
         tmpcart.vz = file.read_f64::<LittleEndian>().unwrap();
         cartesians.push(tmpcart);
    }
        println!("first cartesian {}", cartesians[0].x);
    let mut rads:Vec<f64> = Vec::new();
    let mut tmprad:f64 = 0.0;
    for i in (1..num) {
        tmprad = file.read_f64::<LittleEndian>().unwrap();
        rads.push(tmprad);
        println!("{}", i);
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
        let mut xbin:i32 = ((cartesians[i].x - xmin) * grid.len() as f64 / (xmax - xmin)) as i32;
        let mut ybin:i32 = ((cartesians[i].y - ymin) * grid[0].len() as f64 / (ymax - ymin)) as i32;
        if xbin as usize >= xlen {
            xbin  = xlen as i32 - 1;
        }
        if ybin as usize >= ylen {
            ybin  = ylen as i32 - 1;
        }
        let xylen = grid[xbin as usize][ybin as usize].len();
        grid[xbin as usize][ybin as usize].insert(xylen, i as i32);
   }
   let mut zvec:Vec<f64> = Vec::new();
   for i in (0..cartesians.len()) {
    zvec.push(cartesians[i].z);
   }
   zvec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
   let zmin = zvec[100];
   let zmax = zvec[zvec.len() - 100];
    synthetic_occultation(x, y, theta, phi, cut_theta, scan_length, off_length, beam_size, &binData
                          {xmin: xmin,
                          xmax: xmax,
                          ymin: ymin,
                          ymax: ymax,
                          bin_size: 2e-8,
                          bins: grid,
                          }, zmin, zmax, photon_count, &cartesians, &rads); 
}


fn synthetic_occultation(x: f64, y: f64, theta: f64, phi: f64, cut_theta: f64, scan_length: f64, off_length: f64, beam_size: f64, bin_data: &binData, zmin: f64, zmax: f64, photon_count: i32, cartesians: &Vec<cart>, rad: &Vec<f64>) -> Vec<scan> {
    println!("in synth occ");
    let  r_dir:vec3D = vec3D {
        x: theta.cos() * phi.cos(), 
        y: theta.sin() * phi.cos(),
        z: phi.sin(),
    };
    let dx = cut_theta.cos();
    let dy = cut_theta.sin();
    let mut height = 0 as f64;
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
    let mut sx:f64 = 0.0;
    let mut sy:f64 = 0.0;
    let mut ex:f64 = 0.0;
    let mut ey:f64 = 0.0;
    let mut t:f64 = 0.0;
    let mut rx:f64 = 0.0;
    let mut ry:f64 = 0.0;
    while mx < xend {
        sx = mx;
        sy = y + (cut_theta.tan() * (sx - x));
        ex = sx + scan_length * cut_theta.cos();
        ey = y + cut_theta.tan() * (ex - x);
        let mut photons:Vec<photon> = Vec::new();
        for i in (1..photon_count) {
            let t = rand::random::<f64>();
            let rx = sx + t * (ex - sx) + rand::random::<f64>() * rand::random::<f64>() * beam_size;
            let ry = sy + t * (ey - sy) + rand::random::<f64>() * rand::random::<f64>() * beam_size;
            photons.push(photon {x: rx, y: ry, hit: ray_grid_intersect(ray{r0: vec3D{x: rx, y: ry, z: 0 as f64}, r: r_dir.clone()}, bin_data, zmin, zmax, &cartesians, &rad)})

        }
    ret.push(scan{sx: sx, sy: sy, ex: ex, ey: ey, intensity: (photons.iter_mut().filter(|p|{!p.hit}).count()) as f64 / photon_count as f64, photons:
    photons});
    mx = mx + (scan_length + off_length) * cut_theta.cos();
    }
    return ret;
}

fn ray_grid_intersect(r: ray, bin_data: &binData, zmin: f64, zmax: f64, cartesians: &Vec<cart>, rad: &Vec<f64>) -> bool {
//    println!("in ray_grid_intersect");
    let tmin = (zmin - r.r0.z) / r.r.z;
    let tmax = (zmax - r.r0.z) / r.r.z;
    let xmin = r.r0.x + tmin * r.r.x;
    let xmax = r.r0.x + tmax * r.r.x;
    let ymin = r.r0.y + tmin * r.r.y;
    let ymax = r.r0.y + tmax * r.r.y;
    let minxbin = max(((min(xmin, xmax) - bin_data.xmin) / bin_data.bin_size) - 1 as f64, 0 as f64)
        as i32;
    let maxxbin = min(((max(xmin, xmax) - bin_data.xmin) / bin_data.bin_size) + 1 as f64,
    bin_data.bins.len() as f64 - 1 as f64) as i32;
    let minybin = max(((min(ymin, ymax) - bin_data.ymin) / bin_data.bin_size) - 1 as f64, 0 as f64)
        as i32;
    let maxybin = min(((max(ymin, ymax) - bin_data.ymin) / bin_data.bin_size) + 1 as f64,
    bin_data.bins[0].len() as f64 - 1 as f64) as i32;
    println!("minx: {} miny: {} maxx: {} maxy: {}", minxbin, minybin, maxxbin, maxybin);
    for xbin in (minxbin..maxxbin) {
//      println!("first for");
        for ybin in (maxybin..minybin) {
 //           println!("second for");
            if ray_bin_intersect(r.clone(), xbin, ybin, bin_data, &cartesians, &rad) {
                return true;   
            }
        }
    }
    return false;
}

fn ray_bin_intersect(r: ray, xbintmp: i32, ybintmp: i32, bin_data: &binData, cartesians: &Vec<cart>, rad: &Vec<f64>) -> bool {
//    println!("in ray_bin_intersect");
    let mut xbin = xbintmp;
    let mut ybin = ybintmp;
    while xbin >= bin_data.bins.len() as i32 {
        xbin = xbin - 1;
    }
    while ybin >= bin_data.bins[xbin as usize].len() as i32 {
        ybin = ybin - 1;
    }
    for i in bin_data.bins[xbin as usize][ybin as usize].iter() {
            if (ray_particle_intersect(&r, &cartesians[*i as usize].to_particle(rad[*i as usize]))) {
                return true;
            }
        }
    return false;
}

fn ray_particle_intersect(r: &ray, part: &particle) -> bool {
//    println!("in ray_particle_intersect");
    let p = vec3D{
        x: part.x,
        y: part.y,
        z: part.z,
    };
    let d = vec3D::sub(&r.r0, &p);
    let a = vec3D::dot(&r.r, &r.r);
    let b = vec3D::dot(&r.r, &d) * 2 as f64;
    let c = vec3D::dot(&d, &d) - part.rad * part.rad;
    return b * b - 4 as f64 * a * c >= 0 as f64;
}
/*
fn bin_particles (parts: Vec<particle>) -> binData {
    
    let tmpvec:Vec<Vec<Vec<particle>>> = Vec::new();
    return binData {
    xmin: 0 as f64,
    xmax: 0 as f64,
    ymin: 0 as f64,
    ymax: 0 as f64,
    bin_size: 0 as f64,
    bins: tmpvec,
    }
}
*/
fn min(a: f64, b:f64) -> f64 {
    if a < b {
        return a;
    }
    else {
        return b;
    }
}
fn max(a: f64, b:f64) -> f64 {
    if a > b {
        return a;
    }
    else {
        return b;
    }
}
