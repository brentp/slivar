time ./src/slivarpkg/ddc data/rgp.15.X.bcf data/rgp.ped \
    --chrom "chr15" \
    --exclude /data/human/LCR-hs38.bed.gz \
    --info-fixed "VQSLOD > 0" \
    --info QD \
    --info '^+ClippingRankSum' \
    --info '^+ReadPosRankSum' \
    --info '^+MQRankSum' \
    --info '^+BaseQRankSum' \
    --info '^FS' \
    --info '^ExcessHet' \
    --info DP \
    --info '^SOR' \
    --info VQSLOD \
    --info QUAL \
    --fmt-fixed-het 'AB>=0.2' \
    --fmt-fixed-het 'GQ>=20' \
    --fmt-fixed-hom-ref 'GQ>=20' \
    --fmt-fixed-hom-ref 'AB < 0.015' \

    exit
DONE

time ./src/slivarpkg/ddc ceph/ceph.small.bcf ceph/ceph-subset.ped \
    --chrom "15" \
    --exclude /data/human/LCR-hs37d5.bed.gz \
    --info QD \
    --info '^+ClippingRankSum' \
    --info '^+ReadPosRankSum' \
    --info '^+MQRankSum' \
    --info '^+BaseQRankSum' \
    --info '^FS' \
    --info '^ExcessHet' \
    --info DP \
    --info '^SOR' \
    --info VQSLOD \
    --info QUAL \
    --fmt-fixed-het 'AB>=0.2' \
    --fmt-fixed-het 'GQ>=20' \
    --fmt-fixed-hom-ref 'GQ>=20' \
    --fmt-fixed-hom-ref 'AB < 0.015' \
