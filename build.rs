extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/C/geodesic.c")
        .compile("geographiclib");
}
