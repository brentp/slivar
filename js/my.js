// hi quality variants
function hq(kid, mom, dad) {
    return (kid.DP > 6 && mom.DP > 6 && dad.DP > 6 && mom.GQ > 9 && dad.GQ > 9 && kid.GQ > 9)
}

function hq1(sample) {
  return !sample.unknown && sample.DP > 6 && sample.GQ > 9
}

function hqrv(variant, INFO, af_cutoff) {
  // hi-quality, rare variant.
  return INFO.gnomad_popmax_af < af_cutoff && variant.FILTER == 'PASS'
}

function denovo(kid, mom, dad){
  // check genotypes match expected pattern
  if(!(kid.het && mom.hom_ref && dad.hom_ref && kid.AB > 0.2 && kid.AB < 0.8)){ return false; }
  if(!hq(kid, mom, dad)){ return false; }
  return ((mom.AD[1] + dad.AD[1]) < 2)
}

function x_denovo(kid, mom, dad) {
  if(!(kid.alts >= 1 && mom.hom_ref && dad.hom_ref && kid.AB > 0.3)){ return false; }
  if(!hq(kid, mom, dad)) { return false; }
  if(kid.sex != 'male') { return false; }
  return ((mom.AD[1] + dad.AD[1]) < 2);
}

function hq_ab(sample) {
  // make sure the allele balance is close to 0 for homref, 1 for homalt
  // and intermediate for het
  if(sample.alts == 0){
    return sample.AB < 0.01
  }
  if(sample.alts == 1){
    return sample.AB > 0.2 && sample.AB < 0.8
  }
  if(sample.alts == 2){
    return sample.AB > 0.99
  }
  return false
}

function uniparent_disomy(kid, mom, dad) {
  if(kid.DP < 10 || mom.DP < 10 || dad.DP < 10) { return false; }
  if(kid.GQ < 20 || mom.GQ < 20 || dad.GQ < 20) { return false; }
  if(!(hq_ab(kid) && hq_ab(mom) && hq_ab(dad))){ return false; }
  return (kid.alts == 0 || kid.alts == 2) && ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2));
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
