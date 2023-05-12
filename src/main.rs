#[macro_use] extern crate itertools;

use std::io::{BufWriter, Write};
use std::fs::{File, OpenOptions};
use std::path::PathBuf;
use rug::float::Round;
use rug::Float;
use fraction::GenericFraction;
use core::cmp::Ordering;
use std::collections::BinaryHeap;
use std::cmp::Reverse;
use std::fmt;
use ordered_float::NotNan;
use std::f64::consts::PI;
use rug::float::Constant::Pi as rug_PI;
use itertools::{Itertools, join};
use pico_args;

pub type Fraction = GenericFraction<u64>;
type MinNotNan = Reverse<NotNan<f64>>;

// this is a translation of the existing Python progressive edge search code by Marek and me into Rust
// I will try to optimize little things where possible
// coding in Rust is going kind of slow and challenging
// one big change is variable scope
// also, list comprehensions look less clean when translated literally to Rust

const HELP: &str = "\
Edgecompare

Usage: edgecompare -e MAX_EXP -g MAX_G POLYGON_LIST OUTFILE
  or:  edgecompare -p MAX_P POLYGON_LIST OUTFILE

FLAGS:
  -h, --help            Prints help information

OPTIONS:
  -e, --max-exp MAX_EXP   Sets maximum number of polygons of any given type
  -g, --maxg MAX_G        Sets an optional number
  -p, --maxp MAX_P        Sets maximum number of polygons, counting polygons of all types

ARGS:
  POLYGON_LIST   List of polygon types to include in search
  OUTFILE        Path of output file
";

#[derive(Debug)]
struct AppArgs {
    max_exp: Option<u8>,
    max_g: Option<u8>,
    max_total_exp: Option<u8>,
    polygon_list: (Vec<u8>, Vec<Vec<u8>>),
    outfile: PathBuf,
}


const MAX_IT: usize = 60;
const UNUSED: f64 = -1.1e30;
const PREC: u32 = 169;  // corresponds to Python's mpmath mp.dps = 50
const PRINT_DPS: usize = 50; // show this many significant digits when printing
const M_TOL: f64 = 1e-14;

// somewhat arbitrary thresholds, not sure what is best
const MARGIN_BITS_R: i32 = 4;
const MARGIN_BITS_G: i32 = 10;

fn apply_ratios(m: f64, ratios: &Vec<u8>) -> Vec<f64> {
    let halfedge = m.acosh();
    ratios.into_iter().map(|&r| (r as f64 * halfedge).cosh()).collect()
}

fn apply_ratios_mp(m: &Float, ratios: &Vec<u8>) -> Vec<Float> {
    let halfedge = Float::with_val(PREC, m.acosh_ref());
    ratios.into_iter()
        .map(|r| Float::with_val(PREC, &halfedge * r).cosh())
        .collect()
}

fn vertex_angle(ratios: &Vec<u8>, cosines: &Vec<f64>, ns: &Vec<f64>, m: f64) -> f64 {
    let ms = apply_ratios(m, ratios);
    izip!(cosines, ns, ms)
        .map(|(cs, n, m)| n * 2.0 * (cs/m).asin())
        .sum()
}

fn vertex_angle_mp(ratios: &Vec<u8>, cosines: &Vec<Float>, ns: &Vec<u8>, m: &Float) -> Float {
    let ms = apply_ratios_mp(m, ratios);
    let values: Vec<_> = izip!(cosines, ns, ms)
        .map(|(cs, n, m)| Float::with_val(PREC, cs/m).asin() * 2 * n)
        .collect();
    Float::with_val(PREC, Float::sum(values.iter()))
}

fn angle_fn_mp((ratios, cosines, ns): (&Vec<u8>, &Vec<Float>, &Vec<u8>), m: &Float) -> Float {
    let ms = apply_ratios_mp(m, ratios);
    let values: Vec<_> = izip!(cosines, ns, ms)
        .map(|(cs, n, m)| Float::with_val(PREC, cs/m).asin() * 2 * n)
        .collect();
    Float::with_val(PREC, Float::sum(values.iter())) - 2 * Float::with_val(PREC, rug_PI)
}

fn fmt_edge_length(x: &Float) -> String {
    let (sign, mut s, Some(exp)) = x.to_sign_string_exp_round(10, Some(PRINT_DPS), Round::Nearest) else {
        panic!("Float does not have exponent");
    };
    assert!(!sign);
    assert!(exp <= PRINT_DPS as i32);
    if 0 < exp && exp < PRINT_DPS as i32 {
        s.insert(exp as usize, '.');
    } else if exp == 0 {
        s.insert_str(0, "0.");
    } else {
        s = format!("0.{}{}", "0".repeat((-exp) as usize), s);
    }
    s
}

fn parse_polygon_list(s: &str) -> Result<(Vec<u8>, Vec<Vec<u8>>), &'static str> {
    let mut angles: Vec<(u8, u8)> = if s.chars().any(|x| x.is_ascii_uppercase()) {
        s.split_whitespace().map(|x| parse_angle_symbol_raw(x)).collect::<Result<Vec<_>,_>>()?
    } else {
        s.split("; ")
            .map(|x| parse_ratio_side_range(x))
            .collect::<Result<Vec<_>,_>>()?
            .concat()
    };
    angles.sort_by(|x, y| x.0.cmp(&y.0));

    let mut side_ratios: Vec<u8> = Vec::new();
    let mut polys: Vec<Vec<u8>> = Vec::new();
    for (ratio, group) in &angles.into_iter().group_by(|a| a.0) {
        side_ratios.push(ratio);
        let mut p_list: Vec<u8> = group.map(|(_, p)| p).unique().collect();
        p_list.sort();
        if p_list.contains(&0u8) {
            p_list.retain(|p| *p != 0u8);
            p_list.push(0u8);
        }
        polys.push(p_list);
    }
    Ok((side_ratios, polys))
        // let str_split: Vec<Vec<String>> = s.split("; ")
        //     .map(|s| s.split_whitespace().map(String::from).collect())
        //     .collect();
        // let side_ratios = str_split.iter()
        //     .map(|siter| siter[0].to_string().parse::<u8>())
        //     .collect::<Result<Vec<_>, _>>().map_err(|_| "Unable to parse side ratio")?;
        // let polys = str_split.iter()
        //     .map(|siter| siter[1].split(",")
        //             .map(|x| x.parse::<u8>())
        //             .collect::<Result<Vec<_>, _>>())
        //     .collect::<Result<Vec<_>, _>>().map_err(|_| "Unable to parse polygon side number")?;
        // Ok((side_ratios, polys))
}

fn parse_angle_symbol_raw(s: &str) -> Result<(u8, u8), &'static str> {
    let mut it = s.chars();
    let c = it.next().ok_or("angle symbol parse error")?;
    if c < 'A' || 'Z' < c {
        return Err("Error while parsing angle symbol, expected uppercase letter");
    }
    let ratio: u8 = (c as u32 - 'A' as u32) as u8 + 1u8;

    let sp = it.as_str();
    let p: u8 = if sp == "oo" {
        0u8
    } else {
        sp.parse::<u8>().map_err(|_| "error parsing polygon's number of sides")?
    };

    Ok((ratio, p))
}

fn parse_uint_range(s: &str) -> Result<Vec<u8>, &'static str> {
    let mut result: Vec<u8> = vec![];
    let parts = s.split(",");
    for part in parts {
        if part.contains("..") {
            let part_split = part.split("..").collect::<Vec<_>>();
            if part_split.len() != 2 {
                return Err("Error parsing range, wrong use of ..");
            }
            let a: u8 = part_split[0].parse().map_err(|_| "Error parsing u8")?;
            let b: u8 = part_split[1].parse().map_err(|_| "Error parsing u8")?;
            if a > b {
                return Err("a..b must have a <= b");
            }
            result.extend(a..=b);
        } else {
            result.push(part.parse::<u8>().map_err(|_| "Error parsing u8")?);
        }
    }
    Ok(result)
}

fn parse_ratio_side_range(s: &str) -> Result<Vec<(u8, u8)>, &'static str> {
    let str_split = s.split_whitespace().collect::<Vec<_>>();
    if str_split.len() != 2 {
        return Err("Error while parsing ratio side range, wrong number of spaces");
    }

    let ratios = parse_uint_range(str_split[0])?;
    let nsides = parse_uint_range(str_split[1])?;

    Ok(iproduct!(ratios.into_iter(), nsides.into_iter()).collect::<Vec<(u8, u8)>>())
}

enum SearchLimits {
    GroupLimit { max_exp: u8, max_g: u8 },
    MaxTotalExp(u8),
}

struct ProgressiveEdgeSearch {
    side_ratios: Vec<u8>,
    polys: Vec<Vec<u8>>,
    limits: SearchLimits,
    apro_idx: u8,
    m_tolerance: f64,
    writer: BufWriter<File>,
}

impl ProgressiveEdgeSearch {
    // pub fn new(side_config_str: String, limits: SearchLimits, fname: String) -> Self {
    //     let str_split: Vec<Vec<String>> = side_config_str.split("; ")
    //         .map(|s| s.split_whitespace().map(String::from).collect())
    //         .collect();
    //     let side_ratios = str_split.iter()
    //         .map(|siter| siter[0].to_string().parse::<u8>().unwrap())
    //         .collect::<Vec<_>>();
    //     let polys = str_split.iter()
    //         .map(|siter| siter[1].split(",").map(|x| x.parse::<u8>().unwrap()).collect::<Vec<_>>())
    //         .collect::<Vec<_>>();

    pub fn new(side_ratios: Vec<u8>, polys: Vec<Vec<u8>>, limits: SearchLimits, outfile: PathBuf) -> Self {
        // designate the index used to represent apeirogons
        let apro_idx = polys.iter().flatten().max().unwrap() + 1;

        let file = OpenOptions::new()
            .write(true).create(true)
            .open(outfile)
            .expect("Unable to open file");
        let writer = BufWriter::with_capacity(100, file);

        ProgressiveEdgeSearch { side_ratios, polys, limits, apro_idx, m_tolerance: M_TOL, writer }
    }

    fn is_hyperbolic(&self, atypes: &Vec<AngleType>) -> bool
    {
        let mut invsum = Fraction::from(0);
        let mut n: u64 = 0;
        for a in atypes {
            if a.p < self.apro_idx {
                invsum += Fraction::new(a.mult, a.p);
            }
            n += a.mult as u64;
        }
        return n > 2 && Fraction::new(n-2, 2u8) > invsum;
    }

    fn get_cosines(&self, atypes: &Vec<AngleType>) -> Vec<f64> {
        atypes.into_iter().map(|a| if a.p < self.apro_idx {(PI / (a.p as f64)).cos()} else {1.0}).collect()
    }

    fn get_cosines_mp(&self, atypes: &Vec<AngleType>) -> Vec<Float> {
        atypes.into_iter()
            .map(|a| if a.p < self.apro_idx {
                Float::with_val(PREC, a.p).recip().cos_pi()
            } else {
                Float::with_val(PREC, 1)
            })
            .collect()
    }

    fn get_m(&self, atypes: &Vec<AngleType>) -> f64
    {
        if atypes.len() == 1 {
            self.get_m_single(&atypes[0])
        } else {
            self.get_m_ridder(atypes, 1e-17)
        }
    }

    fn get_edge_precise(&self, atypes: &Vec<AngleType>, guess: Option<&f64>) -> Float {
        if atypes.len() == 1 {
            self.get_edge_precise_single(&atypes[0])
        } else {
            let m = match guess {
                Some(&m) => self.get_m_precise_ridder(atypes, Some(m - 2.0*self.m_tolerance), Some(m + 2.0*self.m_tolerance), false),
                None => self.get_m_precise_ridder(atypes, None, None, false),
            };
            m.acosh() * 2
        }
    }

    fn get_m_single(&self, a: &AngleType) -> f64 {
        let mut m: f64 = (if a.p < self.apro_idx {(PI / (a.p as f64)).cos()} else {1.0}) / (PI / (a.mult as f64)).sin();
        if a.ratio != 1u8 {
            m = (m.acosh()/(a.ratio as f64)).cosh()
        }
        m
    }

    fn get_edge_precise_single(&self, a: &AngleType) -> Float {
        let mut n = if a.p < self.apro_idx {
            Float::with_val(PREC, a.p).recip().cos_pi()
        } else {
            Float::with_val(PREC, 1)
        };
        n /= Float::with_val(PREC, a.mult).recip().sin_pi();
        n = 2 * n.acosh();
        if a.ratio != 1u8 {
            n /= a.ratio;
        }
        n
    }

    fn get_m_ridder(&self, atypes: &Vec<AngleType>, xacc: f64) -> f64 {
        let ratios = atypes.iter().map(|a| a.ratio).collect::<Vec<_>>();
        let cosines = self.get_cosines(atypes);
        let ns = atypes.iter().map(|a| a.mult as f64).collect::<Vec<_>>();
        let rmin = ratios.iter().min().unwrap();
        let rmax = ratios.iter().max().unwrap();
        let sreg = (PI / ns.iter().sum::<f64>()).sin();

        let mut xl = f64::max(1.0, cosines.iter().min_by(|a, b| a.total_cmp(b)).unwrap() / sreg);
        let mut xh = cosines.iter().max_by(|a, b| a.total_cmp(b)).unwrap() / sreg;
        if *rmax != 1 {
            xl = (xl.acosh() / (*rmax as f64)).cosh();
        }
        if *rmin != 1 {
            xh = (xh.acosh() / (*rmin as f64)).cosh();
        }

        // closure
        let angle_c = |x| vertex_angle(&ratios, &cosines, &ns, x) - 2.0*PI;

        let mut fl = angle_c(xl);
        let mut fh = angle_c(xh);
        let mut ans = UNUSED;

        for _ in 0..MAX_IT {
            let xm = (xl + xh)/2.0;
            let fm = angle_c(xm);
            let s = (fm*fm - fl*fh).sqrt();

            if s == 0.0 {
                return ans;
            }
            let xnew = xm + (xm - xl) * (if fl >= fh {1.0} else {-1.0}) * fm / s;
            if (xnew - ans).abs() <= xacc {
                return ans;
            }
            ans = xnew;
            let fnew = angle_c(ans);
            if fnew == 0.0 {
                return ans;
            }
            if fm.copysign(fnew) != fm {
                xl = xm;
                fl = fm;
                xh = ans;
                fh = fnew;
            } else if fl.copysign(fnew) != fl {
                xh = ans;
                fh = fnew;
            } else if fh.copysign(fnew) != fh {
                xl = ans;
                fl = fnew;
            } else {
                unreachable!("Signs of fm, fl, fh, fnew are impossible");
            }
            if (xh-xl).abs() <= xacc {
                return ans;
            }
        }

        println!("Error: Ridders' method exceeded maximum iterations");
        0.0
    }

    fn get_m_precise_ridder(&self, atypes: &Vec<AngleType>, xlo: Option<f64>, xhi: Option<f64>, verbose: bool) -> Float {
        let ratios = atypes.iter().map(|a| a.ratio).collect::<Vec<_>>();
        let cosines = self.get_cosines_mp(atypes);
        let ns = atypes.iter().map(|a| a.mult).collect::<Vec<_>>();

        let (mut xl, mut xh): (Float, Float);
        if let (Some(xl_f), Some(xh_f)) = (xlo, xhi) {
            xl = Float::with_val(PREC, xl_f);
            xh = Float::with_val(PREC, xh_f);
        } else {
            let sreg = Float::with_val(PREC, ns.iter().map(|&x| x as u64).sum::<u64>()).recip().sin_pi();
            xl = match xlo {
                Some(xl_f) => Float::with_val(PREC, xl_f),
                None => {
                    let mut x = Float::with_val(PREC, cosines.iter().min_by(|a, b| a.total_cmp(b)).unwrap()); // min cosine
                    let div = Float::with_val(PREC, &x / &sreg);
                    x = div.max(&Float::with_val(PREC, 1));
                    let rmax = ratios.iter().max().unwrap();
                    if *rmax != 1 {
                        x = (x.acosh() / rmax).cosh();
                    }
                    x
                },
            };
            xh = match xhi {
                Some(xh_f) => Float::with_val(PREC, xh_f),
                None => {
                    let mut x = Float::with_val(PREC, cosines.iter().max_by(|a, b| a.total_cmp(b)).unwrap()); // max cosine
                    x /= sreg;
                    let rmin = ratios.iter().min().unwrap();
                    if *rmin != 1 {
                        x = (x.acosh() / rmin).cosh();
                    }
                    x
                }
            }
        }

        let params = (&ratios, &cosines, &ns);

        let mut fl = angle_fn_mp(params, &xl);
        let mut fh = angle_fn_mp(params, &xh);
        let mut ans = Float::with_val(PREC, UNUSED);
        let (mut xm, mut fm): (Float, Float);

        let xacc: f64 = f64::powi(2.0, -(PREC as i32 - MARGIN_BITS_R));

        if verbose {
            println!("start, f(x) {fl} {fh}");
        }

        for it in 0..MAX_IT {
            xm = Float::with_val(PREC, &xl + &xh) / 2;
            fm = angle_fn_mp(params, &xm);

            // s = (fm*fm - fl*fh).sqrt();
            let mut s = Float::with_val(PREC, fm.square_ref());
            s -= Float::with_val(PREC, &fl * &fh);
            s.sqrt_mut();

            if s.cmp0().unwrap() == Ordering::Equal {
                return ans;
            }

            // xnew = xm + (xm - xl) * (if fl >= fh {1.0} else {-1.0}) * fm / s;
            let mut xnew = if fl >= fh {
                Float::with_val(PREC, &xm - &xl)
            } else {
                Float::with_val(PREC, &xl - &xm)
            };
            xnew *= &fm;
            xnew /= &s;
            xnew += &xm;

            if Float::with_val(PREC, &xnew - &ans).abs() <= xacc {
                return ans;
            }

            ans = xnew;
            let fnew = angle_fn_mp(params, &ans);
            if verbose {
                println!("#{it}, x {ans}, f(xl, xnew, xh) {:+e}, {:+e}, {:+e}", fl.to_f64(), fnew.to_f64(), fh.to_f64());
                println!("   xl {xl} xh {xh}");
            }

            if fnew.cmp0().unwrap() == Ordering::Equal {
                return ans;
            }
            if Float::with_val(PREC, &fm * &fnew) < 0 {
                xl = xm;
                fl = fm;
                xh = Float::with_val(PREC, &ans);
                fh = fnew;
            } else if Float::with_val(PREC, &fl * &fnew) < 0 {
                xh = Float::with_val(PREC, &ans);
                fh = fnew;
            } else if Float::with_val(PREC, &fh * &fnew) < 0 {
                xl = Float::with_val(PREC, &ans);
                fl = fnew;
            } else {
                unreachable!("Signs of fm, fl, fh, fnew are impossible");
            }

            if Float::with_val(PREC, &xh - &xl).abs() <= xacc {
                return ans;
            }
        }

        eprintln!("Error: Ridders' method (multiprecision) exceeded maximum iterations");
        Float::with_val(PREC, 0)
    }

    fn run(&mut self) {
        let mut vertices = BinaryHeap::new();
        for i in 0..self.side_ratios.len() {
            vertices.push(VertexItem::new(
                0.0, 
                vec![IndexedAngleType::new(i as u8, 0u8, 1u8, self)]
            ));
        }

        let mut curmin: Vec<Vec<AngleType>> = vec![];
        let mut curmin_m = 0f64;
        let mut sol = 0u64;
        let report_interval = 1000000u64;
        let t1 = std::time::Instant::now();

        let mut coincidence: u64 = 0;

        while let Some(VertexItem { m: Reverse(n_m), angles: cur_iv }) = vertices.pop() {
            let cur_m: f64 = n_m.into_inner();
            // println!("pop {} {} {}", self.fmt_iatypes(&cur_iv), curmin_m, cur_m);
            
            for new_iv in self.edgebreed(&cur_iv) {
                let new_v = new_iv.iter().map(|x| x.atype).collect();
                // println!("is hyperbolic {}, {:?}", self.is_hyperbolic(&new_v), &new_v);
                vertices.push(
                    VertexItem::new(
                        if self.is_hyperbolic(&new_v) {self.get_m(&new_v)} else {0.0},
                        new_iv
                    )
                )
            }
            sol += 1;
            let cur_v = cur_iv.iter().map(|x| x.atype).collect();
            if sol % report_interval == 0 {
                let t2 = std::time::Instant::now();
                let timeval = t2.duration_since(t1).as_secs_f64();
                println!("{timeval:.2} s, c={sol}, current {}, cosh(e/2)={cur_m}, buffer {}", self.fmt_atypes(&cur_v), vertices.len());
            }
            if cur_m - curmin_m > self.m_tolerance {
                if curmin_m > 0.0 && curmin.len() > 1 {
                    // println!("#{coincidence} vertices: {}", self.format_atypess(&curmin));
                    // coincidence += 1;
                    // if coincidence > 10 {
                    //     break;
                    // }
                    self.processedge(curmin.drain(..).collect(), &curmin_m);
                } else {
                    curmin.clear();
                }
                curmin_m = cur_m;
            }
            curmin.push(cur_v);
        }

        self.writer.flush().expect("Unable to write data");
    }

    fn processedge(&mut self, mut vs: Vec<Vec<AngleType>>, m_min: &f64) {
        // println!("processedge {} {}", self.fmt_atypess(&vs), m_min);

        let edges: Vec<Float> = vs.iter().map(|x| self.get_edge_precise(x, Some(m_min))).collect();
        let e_min = edges.iter().min_by(|a, b| a.total_cmp(b)).unwrap();
        let e_max = edges.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
        let group_threshold: f64 = f64::powi(2.0, -(PREC as i32 - MARGIN_BITS_G));

        if Float::with_val(PREC, e_max-e_min) > group_threshold {
            if vs.len() > 2 {
                println!("? {}", self.fmt_atypess(&vs));
                for edge in edges.iter() {
                    println!("{:?}", edge);
                }

                self.processedge_grouped(vs, edges, group_threshold);
            }
        } else {
            vs.sort();
            self.write_vertices(&vs, &edges[0]);
        }
    }
    
    fn processedge_grouped(&mut self, vs: Vec<Vec<AngleType>>, edges: Vec<Float>, threshold: f64) {
        let mut zipped = vs.into_iter().zip(edges.into_iter())
                                                     .collect::<Vec<(Vec<AngleType>, Float)>>();
        zipped.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        while zipped.len() > 1 {
            let mut i: usize = 1;
            while i < zipped.len() && (Float::with_val(PREC, &zipped[i].1 - &zipped[i-1].1) < threshold) {
                i += 1;
            }
            let (mut v_group, e_group): (Vec<Vec<AngleType>>, Vec<Float>) = zipped.drain(..i).unzip();
            if i > 1 { // want at least two vertices with same edge length
                v_group.sort();
                self.write_vertices(&v_group, &e_group[0]);
            }
        }
    }

    fn write_vertices(&mut self, vs: &Vec<Vec<AngleType>>, edge: &Float) {
        let edge_str = fmt_edge_length(edge);

        println!();
        println!("{}\n{}", &edge_str, self.fmt_atypess(&vs));
        println!();
        
        write!(self.writer, "{}\n{}\n\n", &edge_str, self.fmt_vertex_list(&vs)).expect("Unable to write data");
    }

    fn edgebreed(&self, v: &Vec<IndexedAngleType>) -> Vec<Vec<IndexedAngleType>> {
        let mut children = vec![];
        let &IndexedAngleType{atype: AngleType{ratio: _, p: _, mult: l_mult}, 
            ratio_idx: l_ratioi, p_idx: l_pidx} = v.last().unwrap();

        let advance_pidx = (l_pidx as usize) < self.polys[l_ratioi as usize].len() - 1;

        let (increment_mult, add_next_p, add_later_ratios, replace_p_new_group) = match self.limits {
            SearchLimits::GroupLimit { max_exp, max_g } => {
                let add_group = v.len() < (max_g as usize);

                (l_mult < max_exp,
                l_mult == max_exp && add_group && advance_pidx,
                add_group,
                l_mult > 1 && add_group && advance_pidx)
            },
            SearchLimits::MaxTotalExp(max_total_exp) => {
                let add_polygon = v.iter().map(|x| x.atype.mult).sum::<u8>() < max_total_exp;

                (add_polygon,
                false,
                add_polygon,
                l_mult > 1 && advance_pidx)
            }
        };

        let replace_p = l_mult == 1 && advance_pidx;

        if increment_mult {
            let mut v2 = v.clone();
            // v2.remove(v2.len()-1);
            // v2.push(IndexedAngleType::new(l_ratioi, l_pidx, l_mult+1, self));
            v2.last_mut().unwrap().atype.mult += 1;
            children.push(v2);
        }
        if add_next_p {
            let mut v2 = v.clone();
            v2.push(IndexedAngleType::new(l_ratioi, l_pidx+1, 1, self));
            children.push(v2);
        }
        if add_later_ratios {
            for r in l_ratioi+1..(self.side_ratios.len() as u8) {
                let mut v2 = v.clone();
                v2.push(IndexedAngleType::new(r, 0, 1, self));
                children.push(v2);
            }
        }
        if replace_p {
            let mut v2 = v.clone();
            v2.remove(v2.len()-1);
            v2.push(IndexedAngleType::new(l_ratioi, l_pidx+1, 1, self));
            children.push(v2);
        }
        if replace_p_new_group {
            let mut v2 = v.clone();
            // v2.remove(v2.len()-1);
            // v2.push(IndexedAngleType::new(l_ratioi, l_pidx, l_mult-1, self));
            v2.last_mut().unwrap().atype.mult -= 1;
            v2.push(IndexedAngleType::new(l_ratioi, l_pidx+1, 1, self));
            children.push(v2);
        }

        children
    }

    fn atype(&self, s: &str) -> AngleType {
        let s2: &str;
        let mult: u8;
        if s.contains("^") {
            let parts = s.split("^").collect::<Vec<&str>>();
            s2 = parts[0];
            mult = parts[1].parse::<u8>().unwrap();
        } else {
            s2 = s;
            mult = 1u8;
        }
        let mut it = s2.chars();
        let c = it.next().unwrap();
        assert!('A' <= c);
        assert!(c <= 'Z');
        let ratio: u8 = (c as u32 - 'A' as u32) as u8 + 1u8;

        let sp = it.as_str();
        let p: u8 = if sp == "oo" {
            self.apro_idx
        } else {
            let p2 = sp.parse::<u8>().unwrap();
            if p2 >= self.apro_idx {self.apro_idx} else {p2}
        };

        AngleType { ratio, p, mult }
    }

    fn parse_vertex(&self, s: &str) -> Vec<AngleType> {
        let mut ss = s.trim_start_matches(|c| c == '(' || c == '[');
        ss = ss.trim_end_matches(|c| c == ')' || c == ']');
        let parts: Vec<&str> = if ss.contains(", ") {
            ss.split(", ")
        } else {
            ss.split(",")
        }.collect();
        parts.into_iter().map(|x| self.atype(x)).collect()
    }

    fn fmt_atypes(&self, v: &Vec<AngleType>) -> String {
        format!("({})", join(v.iter().map(|x| format!("{:?}", x.with_source(self))), ","))
    }

    fn fmt_iatypes(&self, v: &Vec<IndexedAngleType>) -> String {
        join(v.iter().map(|x| format!("{:?}", x.with_source(self))), ",")
    }

    fn fmt_atypess(&self, v: &Vec<Vec<AngleType>>) -> String {
        join(v.iter().map(|x| self.fmt_atypes(x)), "; ")
    }

    fn fmt_vertex_list(&self, v: &Vec<Vec<AngleType>>) -> String {
        join(v.iter().map(|x| self.fmt_atypes(x)), "\n")
    }

    fn fmt_iatypess(&self, v: &Vec<Vec<IndexedAngleType>>) -> String {
        join(v.iter().map(|x| self.fmt_iatypes(x)), "; ")
    }
}

/// AngleType instances must be associated with an instance of ProgressiveEdgeSearch, since apeirogons
/// are represented by some "large" integer, typically one larger than the largest finite number of sides 
/// involved in the search
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
struct AngleType {
    ratio: u8,
    p: u8,
    mult: u8,
}

impl Ord for AngleType {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.ratio, self.p, Reverse(self.mult)).cmp(&(other.ratio, other.p, Reverse(other.mult)))
    }
}

impl PartialOrd for AngleType {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

trait Sourceable {
    fn with_source<'a>(&'a self, source: &'a ProgressiveEdgeSearch) -> Sourced<'a, Self> where Self: Sized;
}

struct Sourced<'a, T> where T: Sized {
    val: &'a T,
    source: &'a ProgressiveEdgeSearch,
}

impl<'a> fmt::Debug for Sourced<'a, AngleType> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", std::char::from_u32('A' as u32 + self.val.ratio as u32 - 1).ok_or(fmt::Error)?)?;
        if self.val.p == self.source.apro_idx {
            write!(f, "{}", "oo")?;
        } else {
            write!(f, "{}", self.val.p)?;
        }
        if self.val.mult != 1 {
            write!(f, "^{}", self.val.mult)?;
        }
        Ok(())
    }
}

impl Sourceable for AngleType {
    fn with_source<'a>(&'a self, source: &'a ProgressiveEdgeSearch) -> Sourced<'a, AngleType> {
        Sourced::<'a, AngleType> {
            val: &self,
            source,
        }
    }
}

impl Sourceable for IndexedAngleType {
    fn with_source<'a>(&'a self, source: &'a ProgressiveEdgeSearch) -> Sourced<'a, IndexedAngleType> {
        Sourced::<'a, IndexedAngleType> {
            val: &self,
            source,
        }
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
struct IndexedAngleType {
    atype: AngleType,
    ratio_idx: u8,
    p_idx: u8,
}

impl<'a> fmt::Debug for Sourced<'a, IndexedAngleType> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", std::char::from_u32('A' as u32 + self.val.atype.ratio as u32 - 1).ok_or(fmt::Error)?)?;
        if self.val.atype.p == self.source.apro_idx {
            write!(f, "{}", "oo")?;
        } else {
            write!(f, "{}", self.val.atype.p)?;
        }
        write!(f, "({} {})", self.val.ratio_idx, self.val.p_idx)?;
        if self.val.atype.mult != 1 {
            write!(f, "^{}", self.val.atype.mult)?;
        }
        Ok(())
    }
}

impl IndexedAngleType {
    pub fn new(ratio_idx: u8, p_idx: u8, mult: u8, context: &ProgressiveEdgeSearch) -> Self {
        let p_raw = context.polys[ratio_idx as usize][p_idx as usize];
        IndexedAngleType {
            atype: AngleType{
                ratio: context.side_ratios[ratio_idx as usize],
                p: if p_raw == 0 { context.apro_idx } else { p_raw },
                mult,
            },
            ratio_idx,
            p_idx,
        }
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
struct VertexItem {
    m: MinNotNan,
    angles: Vec<IndexedAngleType>,
}

impl VertexItem {
    pub fn new(mval: f64, angles: Vec<IndexedAngleType>) -> Self {
        VertexItem {
            m: Reverse(NotNan::new(mval).unwrap()),
            angles,
        }
    }
}

fn parse_args() -> Result<AppArgs, pico_args::Error> {
    let mut pargs = pico_args::Arguments::from_env();

    // Help has a higher priority and should be handled separately.
    if pargs.contains(["-h", "--help"]) {
        print!("{}", HELP);
        std::process::exit(0);
    }

    let args = AppArgs {
        max_exp:       pargs.opt_value_from_str(["-e", "--max-exp"])?,
        max_g:         pargs.opt_value_from_str(["-g", "--maxg"])?,
        max_total_exp: pargs.opt_value_from_str(["-p", "--maxp"])?,
        polygon_list:  pargs.free_from_fn(parse_polygon_list)?,
        outfile:       pargs.free_from_os_str(parse_path)?,
    };

    let remaining = pargs.finish();
    if !remaining.is_empty() {
        eprintln!("Warning: unused arguments left: {:?}.", remaining);
    }

    Ok(args)
}

fn parse_path(s: &std::ffi::OsStr) -> Result<std::path::PathBuf, &'static str> {
    Ok(s.into())
}

fn main() {
    let args = match parse_args() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Error: {}.\n", e);
            print!("{}", HELP);
            std::process::exit(1);
        }
    };

    // println!("{:#?}", args);

    let limits: SearchLimits = match (args.max_exp, args.max_g, args.max_total_exp) {
        (Some(_max_exp), Some(_max_g), Some(_max_total_exp)) => {
            eprintln!("Search not implemented for all three limits simultaneously.\n");
            print!("{}", HELP);
            std::process::exit(1);
        },
        (Some(max_exp), Some(max_g), None) => SearchLimits::GroupLimit { max_exp, max_g },
        (None, None, Some(max_total_exp)) => SearchLimits::MaxTotalExp(max_total_exp),
        (_, _, _) => {
            eprintln!("Error: Not enough search limits specified.\n");
            print!("{}", HELP);
            std::process::exit(1);
        }
    };

    let (side_ratios, polys) = args.polygon_list;

    println!("side ratios: {side_ratios:?}");
    println!("polygon # sides: {polys:?}");

    let mut search = ProgressiveEdgeSearch::new(side_ratios, polys, limits, args.outfile);
    search.run();

    // testing, sorry for the mess

//     let side_config_input = "\
// 1 3,0
// 2 0
// 3 3,0
// 4 0
// 6 0
// 10 0
// 12 0";
//     let (max_exp, max_g) = (9, 7);
    // let fname = "test.txt";

    // let side_config_input = "1 3,4,5,8,40,0; 3 0; 5 0; 7 0; 9 0";
    // let limits = SearchLimits::GroupLimit {max_exp: 10, max_g: 6};
    // let fname = "input_p_1_0377.txt";


    // let mut search = ProgressiveEdgeSearch::new(side_config_input.to_string(), 
    //     limits, fname.to_string());
    // // println!("{:?}", search.side_ratios);
    // // println!("{:?}", search.polys);
    // // println!("{}", search.get_m_precise_ridder(&search.parse_vertex("[Joo, B4^3]"), None, None, true));
    // search.run();


    // let mut vertices = BinaryHeap::new();
    // vertices.push(VertexItem::new(3.5, vec![IndexedAngleType::new(0u8, 0u8, 1u8, &search)]));

    // println!("{:?}", vertices);

    // println!("{:?}", search.fmt_iatypess(&search.edgebreed(&vertices.pop().unwrap().angles)));

    // println!("{}", side_config_input);
    // println!("{}", search.is_hyperbolic(&vec![IndexedAngleType::new(2u8, 0u8, 1u8, &search).atype]));

    // println!("{:?}", search.atype("C1^8").with_source(&search));
    // let vert = search.parse_vertex("[A4, B6, Coo^2]");
    // println!("{:?}", vert);
    // println!("{}", search.get_m_ridder(&vert, 1e-17).acosh() * 2.0);
    // println!("{}", search.get_m_single(&search.atype("A3^15")).acosh() * 2.0);
    // println!("{}", search.get_edge_precise(&search.parse_vertex("A3^15"), None));
    // println!("{}", search.get_m_ridder(&search.parse_vertex("[A3, A15^7]"), 1e-17).acosh() * 2.0);
    // println!("{}", search.get_edge_precise(&search.parse_vertex("A3, A15^7"), None));
    // println!("example vertices:");
    // println!("{}", search.get_m_ridder(&search.parse_vertex("[Boo^2, Doo]"), 1e-17));
    // println!("{}", search.get_m_ridder(&search.parse_vertex("[A3^2, B4^3]"), 1e-17));
    // println!("{}", search.get_m_precise_ridder(&search.parse_vertex("[A3^2, B4^3]"), None, None, true));
    // println!("{}", search.get_edge_precise(&search.parse_vertex("[A3^2, B4^3]"), None));

    // test_edge(&search, "(A3^4,Coo,Goo)");
    // test_edge(&search, "(A3,A4^2,A15)");
    // test_edge(&search, "(A4^2,Boo,D4^2,Hoo)");
    // test_edge(&search, "(B3^11,Boo,Coo,Goo^3)");
    // test_edge(&search, "(A3^12,A6^12,B3^6)");

    
    // let w = Float::with_val(PREC, 3u8).recip();


    // // x has a precision of 10 bits
    // let x = Float::with_val(10, 180);
    // // y has a precision of 50 bits
    // let y = Float::with_val(50, Constant::Pi);
    // let incomplete = &x / &y;
    // // z has a precision of 45 bits
    // let z = Float::with_val(45, incomplete);
    // println!("w {:E}, x {}, y {}, z {}", w, x, y, z);
    // assert!(57.295 < z && z < 57.296);

    // test_float(12.3);
    // test_float(1.23);
    // test_float(0.123);
    // test_float(0.0123);

    // let x1: Float = 1.0/Float::with_val(PREC, 3u8);
    // let x2: Float = Float::with_val(PREC, rug_PI);
    // let x3 = Float::with_val(PREC, x2.cosh_ref());
    // println!("x1 {x1}, x2 {x2}, x3 {x3}");
    // let x4 = search.get_edge_precise_single(&search.atype("A5^12"));
    // let x5 = search.get_edge_precise_single(&search.atype("Aoo^4"));
    // println!("x4 {x4}, x5 {x5}");
    
    // let x6 = Float::with_val(15, 3.5);
    // let x7 = Float::with_val(15, &x6 * &2u8);
    // println!("x6 {x6}, x7 {x7}, {}", x7 <= 6.99999);

    // let a: Vec<f64> = vec![2.0, 2.5, -0.5, 1.0, 1.5];
    
    // let maximum = a.iter().max_by(|a, b| a.total_cmp(b)).unwrap() / 1.0;
    // println!("The maximum value was {maximum}.");

}


fn test_edge(search: &ProgressiveEdgeSearch, s: &str) {
    println!("{s}:");
    let mtypes = search.parse_vertex(s);
    let m1 = search.get_m_ridder(&mtypes, 1e-17);
    let m1p = search.get_m_precise_ridder(&mtypes, Some(m1 - search.m_tolerance), Some(m1 + search.m_tolerance), true);
    let edge1 = search.get_edge_precise(&mtypes, Some(&m1));
    println!("{m1}");
    println!("{m1p} {edge1}");
}

fn test_float(x: f64) {
    let my_float = Float::with_val(PREC, x);
    println!("{} {:?}", x, my_float.to_sign_string_exp_round(10, Some(PRINT_DPS), Round::Nearest));
    println!("{}", fmt_edge_length(&my_float));
}