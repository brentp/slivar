// hi quality variants
function hq(kid, mom, dad) {
    return mom.GQ >= 5 && dad.GQ >= 5 && kid.GQ >= 5
}

function hq1(sample) {
  return !sample.unknown && sample.GQ >= 5
}

function denovo(kid, mom, dad){
  // check genotypes match expected pattern
  if(!(kid.het && mom.hom_ref && dad.hom_ref && kid.AB >= 0.2 && kid.AB <= 0.8)){ return false; }
  if(!hq(kid, mom, dad)){ return false; }
  return ((mom.AD[1] + dad.AD[1]) == 0)
}

function x_denovo(kid, mom, dad) {
  if(!(kid.alts >= 1 && mom.hom_ref && dad.hom_ref && kid.AB > 0.3)){ return false; }
  if(!hq(kid, mom, dad)) { return false; }
  if(kid.sex != 'male') { return false; }
  return ((mom.AD[1] + dad.AD[1]) == 0);
}

function recessive(kid, mom, dad) {
  return mom.het && dad.het && kid.alts == 2 && hq(kid, mom, dad)
}

function x_recessive(kid, mom, dad) { 
  return (mom.het && kid.AB > 0.75 && dad.hom_ref && kid.alts >= 1 && hq(kid, mom, dad)
              && kid.sex == 'male' && mom.AB > 0.25 && mom.AB < 0.75)
}

// heterozygous (1 side of compound het)
function solo_ch_het_side(sample) {
  return sample.AB > 0.2 && sample.AB < 0.8 && hq1(sample)
}

function comphet_side(kid, mom, dad) {
  return kid.het && (solo_ch_het_side(mom) != solo_ch_het_side(dad)) && mom.alts != 2 && dad.alts != 2 && solo_ch_het_side(kid) && hq1(mom) && hq1(dad);
}

// assume that mom and kid are affected.
function fake_auto_dom(kid, mom, dad) {
  return kid.het && mom.het && dad.hom_ref && kid.AB >= 0.2 && kid.AB <= 0.8 && mom.AB >= 0.2 && mom.AB <= 0.8 && hq(mom, dad, kid)
}
