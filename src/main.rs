use crate::backend::{Airframe, Atmosphere, TOConf};
use crate::phyeng::RnwCondition;
use crate::planner::{plan, FlapStrategy, PlanInput};
use crate::runway::{Runway, Wind};

mod backend;
mod phyeng;
mod planner;
mod runway;

pub struct Weight {
    pub cg: f32, // % MAC
    pub tow_kg: i32,
}


/*
Vmcg(Minimum Control Speed Ground): Yerde kritik motor arızalandığında aerodinamik dikey kuyruk/rudder ile uçağın pist
merkez hattında tutulabildiği en düşük hız.
Vmca/Vmcw: Havada minimum kontrol hızı.
Vsr(Ref Stall Speed): Seçili konfigürasyondaki min stall hızı.
Vmu(Minimum Unstick Seed): çağın kuyruk sürtmeden emniyetle havalanabileceği asgari geometrik hız.

V2 >= 1.13 * Vsr / 1.20 * Vs
V2 >= 1.10 * Vmca

VR >= 1.05 * Vmca
VR >= 1.05 * Vmu(AEO) / Vmu(OEI)
VR ın alt sınırı 35 feet(agl)deki hızın V2 süratine eşit veya büyük olmasıyla tanımlanır.

V1(ASDA)/(TODA)

ASDA = TORA + SWY
TODA = TORA + CWY

TORA (Take-Off Run Available): Uçağın tekerleklerinin yerle temas ederek hızlanabileceği ve kalkış yapabileceği fiziksel
kaplamalı pist uzunluğudur.

SWY (Stopway - Durma Yolu): Pistin sonunda yer alan, uçağın kalkıştan vazgeçmesi durumunda yapısal hasar görmeden durabilmesi
 için hazırlanmış, yük taşıma kapasitesine sahip ek alandır.

CWY (Clearway - Aşma Yolu): Pistin veya durma yolunun sonunda başlayan, uçağın kalkışın ilk tırmanış safhasında
(ekran yüksekliğine(35 agl) ulaşana kadar) üzerinden emniyetle geçebileceği, engellerden arındırılmış hava/kara sahasıdır.

Bu değerler havalimanı tarafından sağlanmak zorundadır.

Fizik entegrasyonu için ise anlık ivme formülü kullanır.

a = (T - D - u(W-L)) / m

T: Motor İtkisi (Thrust - hıza, sıcaklığa ve basınca bağlı değişir)
D: Aerodinamik Sürükleme (Drag - hızın karesiyle artar)
u: Tekerlek sürtünme katsayısı (Kuru, ıslak veya karlı zemin durumuna göre değişir)
W: Uçağın Ağırlığı (Weight)
L: Taşıma Kuvveti (Lift - hızlandıkça artar, tekerleklere binen yükü azaltır)
m: Kütle

Bu formül anlık integre edilerek katedilen mesafe bulunur( s = intg(v/a dv))

Not:Uçak bu verileri anlık olarak hesaplayamiyacağına göre bu verileri simüle edip hesaplamamamız gerekir.
Bunun için v1 değerini > vmcg olacak şekilde en küçük şekilde seçip interolasyon ile itere ede ede aşşağıda belirilen
ASDR Ve TODR değerlerine uygun hale getirmemiz gerekir böylece v1 değerimizi hesaplarız.

ASDR
Aşşağıdaki kat edilen mesafelerin toplamlarıdır.
1) Fren bırakılan mesafe ile Vef(Engine Failure) arasındaki süre.
2) Reaksiyon süresi dahilinde kat edilen mesafe.(Burda CS-25 standartalrına göre reaksiyon süresi 1 saniye belirlenmiş ancak güvenlik marjı olması için 2 saniye koyuyoruz)
3) Yavaşlama süreci içinde kat edilen mesafe.

TODR
Burda iki farklı hesap vardır ikisi de ayrı ayrı hesaplanır hangisi uzunsa o kabul edilir.
1) OEI durumudur Vef hızında motor arızalanır Vr geçer rotasyon yapar bu rotasyon sürecinde screen height(35 agl) a ulaşana kadarki kat ettiği mesafedir.
2) AEO durumudur uçak hızlanır yine screen height a kadarki mesafesi alınır. Ancak güvenlik marjı için 1.15 ile çarpılır.

Burdaki hesap ASDR ile TODR ı olabildiğince yakınlaştırmakdır v1 değeri böyle bulunur.
Ek olarak TODR <= TODA ve ASDR <= ASDA olmak zorundadır.

Vef >= Vmcg

V1 ile Vef farkı aradaki güvenlik tamponudur Vef motorun bozulduğu hız V1 ise eylem alındığı hızdır.
Dolayısıyla V1 = Vef + (tek motorla belirtilen saniyede(en az 1 olcak şekilde) alınan hız) dır.

Burda hesaplama için gerekli adımlar

1) Vef = Vmcg olarak başlanır
2) Belirtilen formülden(99) ilk v1 hızı çekilir.
3) Bu V1 değeri ile TODR ve ASDR hesaplanır.
4) Eğer ASDR < TODR ise Vef arttırılır(Çünkü pist marjımız vardır) değilse ve TODR < ASDR ise Vef düşürülür(imkansıza yakın ihtimal)
ve eğer ASDR = TODR ise Vef değeri bulunur dolayısıyla V1 de bulunmuş olur.


 */



fn main() {
    let w = Weight { cg: 25.0, tow_kg: 65819 };
    let toconf = TOConf { Packs: false, EngAI: false, WingAI: false };
    let (qnh_hpa, oat_c) = (1013.0, 15.0);

    let rwy = Runway {
        heading_deg: 0.0,
        tora_m: 3000.0,
        swy_m: 0.0,
        cwy_m: 0.0,
        elevation_ft: 0.0,
        slope_pct: 0.0,
    };
    let wind = Wind { dir_deg: 180.0, speed_kt: 21.0 };
    let limits = rwy.limits(&wind);

    let atm = match Atmosphere::new(qnh_hpa, oat_c, rwy.elevation_ft) {
        Ok(a) => a,
        Err(e) => {
            println!("{e}");
            return;
        }
    };

    let path = std::env::args().nth(1).unwrap_or_else(|| "--path to json ".to_string());
    let airframe = match Airframe::get_from_json(&path) {
        Ok(a) => a,
        Err(e) => {
            println!("{e}");
            return;
        }
    };

    let input = PlanInput {
        airframe: &airframe,
        atmosphere: &atm,
        weight: &w,
        toconf: &toconf,
        rnw: RnwCondition::Dry,
        runway: &limits,
        strategy: FlapStrategy::MaxFlex,
        allow_flex: true,
    };

    match plan(&input) {
        Ok(r) => {
            let b = &r.best;
            let t = &b.takeoff;
            let flex = b.flex_temp_c.map_or("TOGA".to_string(), |f| format!("{f} °C"));
            let trim = b.pitch_trim.map_or("-".to_string(), |p| p.to_string());
            println!("FLAP {}  FLEX {}  TRIM {}", b.flap, flex, trim);
            println!("V1 {:.0}  VR {:.0}  V2 {:.0}", t.v1_kt, t.vr_kt, t.v2_kt);
            println!("TODR {:.0} m  ASDR {:.0} m", t.todr_m, t.asdr_m);
            println!("GRADIENT %{:.2}", b.oei_gradient_pct);
        }
        Err(e) => println!("ERROR: {e}"),
    }
}
