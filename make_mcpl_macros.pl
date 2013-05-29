#!/usr/bin/perl -w


open(RUN_LIST,"run_list.dat");
@runs = <RUN_LIST>;

$outfile = "sn.mac";

open(OUT, ">$outfile");

printf OUT "/glg4debug/glg4param omit_muon_processes  1.0 \n";
printf OUT "/glg4debug/glg4param omit_hadronic_processes  1.0 \n";
printf OUT "\n";
printf OUT "/rat/db/set DETECTOR geo_file \"geo/snoplus.geo\" \n";
printf OUT "\n";
printf OUT "/rat/db/set DAQ nhit_thresh 10.0 \n";
printf OUT "/run/initialize \n";
printf OUT "/rat/proc frontend \n";
printf OUT "/rat/proc trigger \n";
printf OUT "/rat/proc eventbuilder \n";
printf OUT "/rat/proc count \n";
printf OUT "/rat/procset update 1 \n";
printf OUT "/rat/proc fitcentroid \n";
printf OUT "/rat/proc outroot \n";
printf OUT "/rat/procset file \"sn_out.root\" \n";
printf OUT "\n";
printf OUT "\n";

foreach $run(@runs) {

    print "working on run ",$run, "\n";

    open(MCPL,"<$run");
    
    $line = 0;
    while(<MCPL>){
	if($line > 32)
	{
	    ($num, $id, $x, $y, $z, $t, $e, $px, $py, $pz) = split;
	    $num = $num;  #just so you dont get a warning ;
	    if ($id == 20){
		$id2 = "e-";
	    }
	    if ($id == 81){
		$id2 = "neutron";
	    }
	    if ($id == 80){
		$id2 = "proton";
	    }
	    if ($id == 2){
		$id2 = "gamma";
	    }
	   
	    printf OUT "/generator/add combo gun3:point \n";
	    printf OUT "/generator/pos/set $x $y $z \n";
	    printf OUT "/generator/vtx/set $id2 $px $py $pz $e $t \n";
	    printf OUT "/generator/rate/set 1 \n";
	    printf OUT "/run/beamOn \n";
	    
	    printf OUT "\n";
	    printf OUT "\n";
	}
	$line++;
    }
    printf OUT "exit\n";
}

