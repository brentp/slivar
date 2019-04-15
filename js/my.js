// hi quality variants
function hq(kid, mom, dad) {
    return (kid.DP > 7 && mom.DP > 7 && dad.DP > 7 
         && mom.GQ > 10 && dad.GQ > 10 && kid.GQ > 10)
}

function hqrv(variant, INFO, af_cutoff) {
	// hi-quality, rare variant.
	return (!('gnomad_popmax_af_filter' in INFO) && variant.FILTER == 'PASS'
         && (!variant.is_multiallelic) && INFO.gnomad_popmax_af_controls < af_cutoff)
}
function denovo(kid, mom, dad){
  if(!hq(kid, mom, dad)){ return false; }
  // check genotypes match expected pattern
  if(!(kid.alts == 1 && mom.alts == 0 && dad.alts == 0)){ return false; }
  return ((mom.AD[1] + dad.AD[1]) < 2) && kid.AB > 0.2 && kid.AB < 0.8;
}

function x_denovo(kid, mom, dad) {
  if(!hq(kid, mom, dad)) { return false; }
  if(kid.sex != 'male') { return false; }
  if(!(mom.alts == 0 && dad.alts == 0 && kid.alts >= 1 && kid.AB > 0.3)){ return false; }
  return ((mom.AD[1] + dad.AD[1]) < 2);
}

function hq_ab(sample) {
  // make sure the allele balance is close to 0 for homref, 1 for homalt
  // and intermediate for het
  if(sample.alts == 0){
    return sample.AB < 0.01
  }
  if(sample.alts == 1){
    return sample.AB > 0.3 && sample.AB < 0.7
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
  return hq(kid, mom, dad) && mom.alts == 1 && dad.alts == 1 && kid.alts == 2
}

function x_recessive(kid, mom, dad) {
  return (hq(kid, mom, dad) && mom.alts == 1 && dad.alts == 0 && kid.alts >= 1
              && kid.sex == 'male' && mom.AB > 0.25 && mom.AB < 0.75 && kid.AB > 0.75)
}

function lenient_denovo(kid, mom, dad){
  if(!hq(kid, mom, dad)){ return false; }
  // check genotypes match expected pattern
  if(!(kid.alts == 1 && mom.alts == 0 && dad.alts == 0)){ return false; }
  return ((mom.AD[1] + dad.AD[1]) < 2) && kid.AB > 0.2 && kid.AB < 0.8;
}

function lenient_ar(kid, mom, dad) {
	 return hq(kid, mom, dad) &&
       (mom.alts == 1 || dad.alts == 1) && kid.alts >= 1 && mom.alts != -1 && dad.alts != -1 &&
        dad.alts != 2 && mom.alts != 2
        && (kid.AB > 0.25 && kid.AB < 0.75)
}
