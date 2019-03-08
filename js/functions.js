function trio_denovo(kid, dad, mom) {
    // alts are 0:hom_ref, 1:het, 2:hom_alt, -1:unknown
    if(!(kid.alts == 1 && mom.alts == 0 && dad.alts == 0)){ return false}
    // sufficient depth in all
    if(kid.DP < 7 || mom.DP < 7 || dad.DP < 7) { return false; }
    // no evidence for alternate in the parents
    if((mom.AD[1] + dad.AD[1]) > 0) { return false; }
    // check the kid's allele balance.
    if(kid.AB < 0.2 || kid.AB > 0.8) { return false; }
    return true
}

function trio_autosomal_dominant(kid, dad, mom) {
    // affected samples must be het.
    if(!(kid.affected == (kid.alts == 1) && mom.affected == (mom.alts == 1) && dad.affected == (dad.alts == 1))) { return false; }
    if(kid.DP < 7 || mom.DP < 7 || dad.DP < 7) { return false; }
    return kid.affected
}

function trio_autosomal_recessive(kid, dad, mom) {
    if(!(kid.affected == (kid.alts == 2) && mom.affected == (mom.alts == 2) && dad.affected == (dad.alts == 2))) { return false; }
    if(kid.DP < 7 || mom.DP < 7 || dad.DP < 7) { return false; }
    if(kid.GQ < 10 || mom.GQ < 10 || dad.GQ < 10) { return false; }
    return kid.affected
}
