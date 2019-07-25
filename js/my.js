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

// functions to be use with --family-expr

function segregating_dominant_x(s) {
  // this is an internal function only called after checking sample quality and on X
  if(!s.affected) { return s.hom_ref }
  if(s.sex == "male") {
    for(var i=0; i < s.kids.length; s++) {
      var kid = s.kids[i];
      // kids of affected dad must be affected.
      if(!kid.affected) { return false; }
    }
    // mom of affected male must be affected.
    if(s.mom && !(s.mom.affected && s.mom.het)){ return false; }

    return s.hom_alt || (s.het && s.AB < 0.2 && s.AB > 0.8)
  }
  if(s.sex != "female"){return false; }
  // this block enforces inherited dominant, but not find de novos
  if(s.mom || s.dad) {
    if(!((s.mom && s.mom.affected && s.mom.het) || (s.dad && s.dad.affected))) { return false;}
  }
  return s.het && s.AB > 0.2 && s.AB < 0.8;
}


function hom_ref_parent(s) {
	return s.dad && s.dad.hom_ref && s.mom && s.mom.hom_ref
}

function segregating_dominant(s) {
  if(s.GQ < 10 || s.unknown) { return false; }
  if(variant.CHROM == "chrX" || variant.CHROM == "X") { return segregating_dominant_x(s); }
  if (s.affected) {
     return s.het && s.AB >= 0.2 && s.AB <= 0.8
  }
  return s.hom_ref && s.AB < 0.01
}


function segregating_recessive_x(s) {
  // this is an internal function only called after checking sample quality and on X
  if(s.sex == "female") {
    return s.affected == s.hom_alt;
  } else if (s.sex == "male") {
    if (s.affected && s.het && hom_ref_parent(s)) { return false; }
    return s.affected == (s.het || s.hom_alt);
  } else {
    return false;
  }
}

function segregating_recessive(s) {
  if(s.GQ < 10 || s.unknown) { return false; }
  if(variant.CHROM == "chrX" || variant.CHROM == "X") { return segregating_recessive_x(s); }
  if(s.affected){
    return s.hom_alt
  }
  return s.het || s.hom_ref
}

// this function is used internally. called from segregating de novo.
// we already know it's on chrX
function segregating_denovo_x(s) {
  if(s.sex == "female") {
    if(s.affected) { return s.het && 0.2 <= s.AB && s.AB <= 0.8 && hom_ref_parent(s) }
    return s.hom_ref
  }
  if(s.sex == "male") {
    if(s.affected)  { return (s.het || s.hom_alt) && hom_ref_parent(s) }	
    return s.hom_ref
  }
  return false
}

function segregating_denovo(s) {
  if (s.GQ < 10 || s.unknown) { return false; }
  if (variant.CHROM == "chrX" || variant.CHROM == "X") { return segregating_denovo_x(s); }
  if (!s.affected) { return s.hom_ref }
  if (s.hom_alt) { return false }
  return s.het && s.AB >= 0.2 && s.AB <= 0.8
}
