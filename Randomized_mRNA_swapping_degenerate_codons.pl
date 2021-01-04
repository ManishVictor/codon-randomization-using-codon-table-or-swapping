#!/usr/bin/perl -w
open (f1,"CDS.fasta");#Input coding sequence file (CDS). The file shoud be in fasta format with only the gene name. eg >YAL007:Mitochondrial gene regulating the water level(Not allowed).>YAL007(Allowed) 
@li=<f1>;
close(f1);
$l=@li;
open (f1,">randomised_sequence.fasta");#The name of the output file having the sequences with randomised synonymous codons within the gene. 

chomp ($li[0]);$na=$li[0];$se='';
$gga=0;$ggt=0;$ggg=0;$ggc=0;$gaa=0;$gat=0;$gag=0;$gac=0;$gca=0;$gct=0;$gcg=0;$gcc=0;$gta=0;$gtt=0;$gtg=0;$gtc=0;
$aga=0;$agt=0;$agg=0;$agc=0;$aaa=0;$aat=0;$aag=0;$aac=0;$aca=0;$act=0;$acg=0;$acc=0;$ata=0;$att=0;$atg=0;$atc=0;
$tga=0;$tgt=0;$tgg=0;$tgc=0;$taa=0;$tat=0;$tag=0;$tac=0;$tca=0;$tct=0;$tcg=0;$tcc=0;$tta=0;$ttt=0;$ttg=0;$ttc=0;
$cga=0;$cgt=0;$cgg=0;$cgc=0;$caa=0;$cat=0;$cag=0;$cac=0;$cca=0;$cct=0;$ccg=0;$ccc=0;$cta=0;$ctt=0;$ctg=0;$ctc=0;$ta=0;
for ($i=1;$i<$l;++$i)
   {chomp ($li[$i]);if ($li[$i]=~/>/){

   $tl=length($se);@chars = split(//, $se);for ($j=0;$j<$tl;) {$cod='';$cod=$chars[$j];++$j;$cod.=$chars[$j];++$j;$cod.=$chars[$j];++$j;++$ta;
	if ($cod eq 'GCA') {$aa.='A';$Aseq.=$cod;++$A;} elsif ($cod eq 'GCC') {$aa.='A';$Aseq.=$cod;++$A;}elsif ($cod eq 'GCG') {$aa.='A';$Aseq.=$cod;++$A;}elsif ($cod eq 'GCT') {$aa.='A';$Aseq.=$cod;++$A;} 
elsif ($cod eq 'GTA') {$aa.='V';$Vseq.=$cod;++$A;} elsif ($cod eq 'GTC') {$aa.='V';$Vseq.=$cod;++$A;} elsif ($cod eq 'GTG') {$aa.='V';$Vseq.=$cod;++$A;} elsif ($cod eq 'GTT') {$aa.='V';$Vseq.=$cod;++$A;}             
elsif ($cod eq 'GGA') {$aa.='G';$Gseq.=$cod;++$G;} elsif ($cod eq 'GGC') {$aa.='G';$Gseq.=$cod;++$G;}elsif ($cod eq 'GGG') {$aa.='G';$Gseq.=$cod;++$G;}elsif ($cod eq 'GGT') {$aa.='G';$Gseq.=$cod;++$G;}
elsif ($cod eq 'GAA') {$aa.='E';$Eseq.=$cod;++$E;} elsif ($cod eq 'GAC') {$aa.='D';$Dseq.=$cod;++$D;}elsif ($cod eq 'GAG') {$aa.='E';$Eseq.=$cod;++$E;}elsif ($cod eq 'GAT') {$aa.='D';$Dseq.=$cod;++$D;}
elsif ($cod eq 'CCA') {$aa.='P';$Pseq.=$cod;++$P;} elsif ($cod eq 'CCC') {$aa.='P';$Pseq.=$cod;++$P;}elsif ($cod eq 'CCG') {$aa.='P';$Pseq.=$cod;++$P;}elsif ($cod eq 'CCT') {$aa.='P';$Pseq.=$cod;++$P;}
elsif ($cod eq 'CTA') {$aa.='L';$Lseq.=$cod;++$L;} elsif ($cod eq 'CTC') {$aa.='L';$Lseq.=$cod;++$L;}elsif ($cod eq 'CTG') {$aa.='L';$Lseq.=$cod;++$L;}elsif ($cod eq 'CTT') {$aa.='L';$Lseq.=$cod;++$L;}
elsif ($cod eq 'CGA') {$aa.='R';$Rseq.=$cod;++$R;} elsif ($cod eq 'CGC') {$aa.='R';$Rseq.=$cod;++$R;}elsif ($cod eq 'CGG') {$aa.='R';$Rseq.=$cod;++$R;}elsif ($cod eq 'CGT') {$aa.='R';$Rseq.=$cod;++$R;}
elsif ($cod eq 'CAA') {$aa.='Q';$Qseq.=$cod;++$Q;} elsif ($cod eq 'CAC') {$aa.='H';$Hseq.=$cod;++$H;}elsif ($cod eq 'CAG') {$aa.='Q';$Qseq.=$cod;++$Q;}elsif ($cod eq 'CAT') {$aa.='H';$Hseq.=$cod;++$H;}
elsif ($cod eq 'TCA') {$aa.='S';$Sseq.=$cod;++$S;} elsif ($cod eq 'TCC') {$aa.='S';$Sseq.=$cod;++$S;}elsif ($cod eq 'TCG') {$aa.='S';$Sseq.=$cod;++$S;}elsif ($cod eq 'TCT') {$aa.='S';$Sseq.=$cod;++$S;}
elsif ($cod eq 'TTA') {$aa.='L';$Lseq.=$cod;++$L;} elsif ($cod eq 'TTC') {$aa.='F';$Fseq.=$cod;++$F;}elsif ($cod eq 'TTG') {$aa.='L';$Lseq.=$cod;++$L;}elsif ($cod eq 'TTT') {$aa.='F';$Fseq.=$cod;++$F;}
elsif ($cod eq 'TGC') {$aa.='C';$Cseq.=$cod;++$C;}elsif ($cod eq 'TGG') {$aa.='W';$Wseq.=$cod;++$W;}elsif ($cod eq 'TGT') {$aa.='C';$Cseq.=$cod;++$C;}
elsif ($cod eq 'TAC') {$aa.='Y';$Yseq.=$cod;++$Y;}elsif ($cod eq 'TAT') {$aa.='Y';$Yseq.=$cod;++$Y;}
elsif ($cod eq 'ACA') {$aa.='T';$Tseq.=$cod;++$T;} elsif ($cod eq 'ACC') {$aa.='T';$Tseq.=$cod;++$T;}elsif ($cod eq 'ACG') {$aa.='T';$Tseq.=$cod;++$T;}elsif ($cod eq 'ACT') {$aa.='T';$Tseq.=$cod;++$T;}
elsif ($cod eq 'ATA') {$aa.='I';$Iseq.=$cod;++$I;} elsif ($cod eq 'ATC') {$aa.='I';$Iseq.=$cod;++$I;}elsif ($cod eq 'ATG') {$aa.='M';$Mseq.=$cod;++$M;}elsif ($cod eq 'ATT') {$aa.='I';$Iseq.=$cod;++$I;}
elsif ($cod eq 'AGA') {$aa.='R';$Rseq.=$cod;++$R;} elsif ($cod eq 'AGC') {$aa.='S';$Sseq.=$cod;++$S;}elsif ($cod eq 'AGG') {$aa.='R';$Rseq.=$cod;++$R;}elsif ($cod eq 'AGT') {$aa.='S';$Sseq.=$cod;++$S;}
elsif ($cod eq 'AAA') {$aa.='K';$Kseq.=$cod;++$K;} elsif ($cod eq 'AAC') {$aa.='N';$Nseq.=$cod;++$N;}elsif ($cod eq 'AAG') {$aa.='K';$Kseq.=$cod;++$K;}elsif ($cod eq 'AAT') {$aa.='N';$Nseq.=$cod;++$N;}
else {$stop=$cod;}
}

$ala='';@new=();my @old = $Aseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$alas.=$new[$j];}@As= split(/(.{1,3})/, $alas);$al=1;
$val='';@new=();my @old = $Vseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$vals.=$new[$j];}@Vs= split (/(.{1,3})/, $vals);$va=1;
$gly='';@new=();my @old = $Gseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$glys.=$new[$j];}@Gs= split(/(.{1,3})/, $glys);$gy=1;
$asp='';@new=();my @old = $Dseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$asps.=$new[$j];}@Ds= split(/(.{1,3})/, $asps);$ap=1;
$pro='';@new=();my @old = $Pseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$pros.=$new[$j];}@Ps= split(/(.{1,3})/, $pros);$prx=1;
$leu='';@new=();my @old = $Lseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$leus.=$new[$j];}@Ls= split(/(.{1,3})/, $leus);$le=1;
$arg='';@new=();my @old = $Rseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$args.=$new[$j];}@Rs= split(/(.{1,3})/, $args);$ar=1;
$glu='';@new=();my @old = $Eseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$glus.=$new[$j];}@Es= split(/(.{1,3})/, $glus);$gu=1;
$gln='';@new=();my @old = $Qseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$glns.=$new[$j];}@Qs= split(/(.{1,3})/, $glns);$gn=1;
$cys='';@new=();my @old = $Cseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$cyss.=$new[$j];}@Cs= split(/(.{1,3})/, $cyss);$cy=1;
$phe='';@new=();my @old = $Fseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$phes.=$new[$j];}@Fs= split (/(.{1,3})/, $phes);$ph=1;
$his='';@new=();my @old = $Hseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$hiss.=$new[$j];}@Hs= split(/(.{1,3})/, $hiss);$hi=1;
$ile='';@new=();my @old = $Iseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$iles.=$new[$j];}@Is= split (/(.{1,3})/, $iles);$il=1;
$lys='';@new=();my @old = $Kseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$lyss.=$new[$j];}@Ks= split(/(.{1,3})/, $lyss);$ly=1;
$met='';@new=();my @old = $Mseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$mets.=$new[$j];}@Ms= split (/(.{1,3})/, $mets);$me=1;
$asn='';@new=();my @old = $Nseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$asns.=$new[$j];}@Ns= split(/(.{1,3})/, $asns);$as=1;
$ser='';@new=();my @old = $Sseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$sers.=$new[$j];}@Ss= split (/(.{1,3})/, $sers);$sl=1;
$thr='';@new=();my @old = $Tseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$thrs.=$new[$j];}@Ts= split(/(.{1,3})/, $thrs);$th=1;
$trp='';@new=();my @old = $Wseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$trps.=$new[$j];}@Ws= split (/(.{1,3})/, $trps);$tr=1;
$tyr='';@new=();my @old = $Yseq =~ /(.{1,3})/g; while (@old) {push(@new, splice(@old, rand @old, 1));}$w=@new;for ($j=0;$j<$w;++$j) {$tyrs.=$new[$j];}@Ys= split (/(.{1,3})/, $tyrs);$ty=1;



@am=split (//,$aa);$w=@am;$ns='';
for ($k=0;$k<$w;++$k)
{$pr=$am[$k];$cod='';

if ($pr eq 'A') {$cod=$As[$al];$al=$al+2;}elsif ($pr eq 'V') {$cod=$Vs[$va];$va=$va+2;}elsif ($pr eq 'G') {$cod=$Gs[$gy];$gy=$gy+2;}
elsif ($pr eq 'D') {$cod=$Ds[$ap];$ap=$ap+2;}elsif ($pr eq 'P') {$cod=$Ps[$prx];$prx=$prx+2;}elsif ($pr eq 'R') {$cod=$Rs[$ar];$ar=$ar+2;}
elsif ($pr eq 'E') {$cod=$Es[$gu];$gu=$gu+2;}elsif ($pr eq 'C') {$cod=$Cs[$cy];$cy=$cy+2;}elsif ($pr eq 'F') {$cod=$Fs[$ph];$ph=$ph+2;}
elsif ($pr eq 'H') {$cod=$Hs[$hi];$hi=$hi+2;}elsif ($pr eq 'I') {$cod=$Is[$il];$il=$il+2;}elsif ($pr eq 'K') {$cod=$Ks[$ly];$ly=$ly+2;}
elsif ($pr eq 'M') {$cod=$Ms[$me];$me=$me+2;}elsif ($pr eq 'N') {$cod=$Ns[$as];$as=$as+2;}elsif ($pr eq 'S') {$cod=$Ss[$sl];$sl=$sl+2;}
elsif ($pr eq 'T') {$cod=$Ts[$th];$th=$th+2;}elsif ($pr eq 'W') {$cod=$Ws[$tr];$tr=$tr+2;}elsif ($pr eq 'Y') {$cod=$Ys[$ty];$ty=$ty+2;}
elsif ($pr eq 'L') {$cod=$Ls[$le];$le=$le+2;}elsif ($pr eq 'Q') {$cod=$Qs[$gn];$gn=$gn+2;}
$ns.=$cod;

}

$ns.=$stop;
print f1 "$na\n$ns\n";
$na=$li[$i];$se='';$ns='';$aa='';$w='';$cod='';$Aseq='';$Cseq='';$Dseq='';$Eseq='';$Fseq='';$Gseq='';$Hseq='';$Iseq='';$Kseq='';$Lseq='';$Mseq='';$Nseq='';$Pseq='';$Qseq='';$Rseq='';$Sseq='';$Tseq='';$Vseq='';$Wseq='';$Yseq='';$ala='';$val='';$gly='';$asp='';$pro='';$arg='';$glu='';$cys='';$phe='';$his='';$ile='';$lys='';$met='';$asn='';$ser='';$thr='';$trp='';$tyr='';$leu='';$gln='';$al='';$va='';$gy='';$ap='';$prx='';$ar='';$gu='';$cy='';$ph='';$hi='';$il='';$ly='';$me='';$as='';$se='';$th='';$tr='';$ty='';$le='';$gn='';$pr='';$A='';$C='';$D='';$E='';$F='';$G='';$H='';$I='';$K='';$L='';$M='';$N='';$P='';$Q='';$R='';$S='';$T='';$V='';$W='';$Y='';$alas='';$vals='';$glys='';$asps='';$pros='';$args='';$glus='';$cyss='';$phes='';$hiss='';$iles='';$lyss='';$mets='';$asns='';$sers='';$thrs='';$trps='';$tyrs='';$leus='';$glns='';
}
else {$se.=$li[$i];}
}

