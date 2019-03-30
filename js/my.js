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

function recessive(kid, mom, dad) {
  return hq(kid, mom, dad) && mom.alts == 1 && dad.alts == 1 && kid.alts == 2
}

function x_recessive(kid, mom, dad) {
  return (hq(kid, mom, dad) && mom.alts == 1 && dad.alts == 0 && kid.alts == 1
              && kid.sex == 'male' && mom.AB > 0.2 && mom.AB < 0.8 && kid.AB > 0.3)
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
