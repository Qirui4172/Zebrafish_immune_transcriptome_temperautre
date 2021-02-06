#!/usr/bin/perl -w

use strict;
use warnings;

#=============================================================================================
sub Usage{
	print "Usage: perl overlapTargets.pl [targetList1] [targetList2] [targetList3] [shortName1] [shortName2] [shortName3] [minOverlap] [outputName]

	Parameters:
	[targetList1]	Target gene list1, only using the first two columns: dre-miR-122-5p	ENSDARG00000000442	slc39a13
	[targetList2]	Target gene list2, same as above
	[targetList3]	Target gene list3, same as above
	[shortName1]	Short name for targetlist1 file
	[shortName2]	Short name for targetlist2 file
	[shortName3]	Short name for targetlist3 file
	[minOverlap]	Minimal overlapping number of target lists, 1, 2 or 3. If set to 1 target genes appeared in either of three target lists will be isolated, if set to 2 targets overlapped in >=2 lists will be isolated, and so on.
	[outputName]	Output file name of overlapped target genes

	Function:
	Isolate overlapped target genes of the same miRNA from three targetList files.

	Example:
	perl overlapTargets.pl mirnaTargets_hybrid.txt mirnaTargets_miranda.txt mirnaTargets_targetscan.txt hybrid miranda targetscan 2 mirnaTargets_overlap2.txt

	Contact: Qirui Zhang (qirui.zhang\@med.lu.se)
	Date: 05-01-2021\n";
}

if(@ARGV!=8){Usage();exit;}


#=============================================================================================
my ($list1, $list2, $list3, $name1, $name2, $name3, $minoverlap, $output)=@ARGV;
my %targetList;

open (LIST1, "$list1") or die "$!\n";
while(<LIST1>){
	chomp;
	my @tmp=(split/\t/, $_); # dre-miR-122-5p ENSDARG00000000442  slc39a13
	my $key=$tmp[0]."_".$tmp[1];
	my $value=$_."\t".$name1."\t1";
	$targetList{$key}=$value;
}
close LIST1;


open (LIST2, "$list2") or die "$!\n";
while (<LIST2>){
	chomp;
	my @tmp=(split/\t/, $_); # dre-miR-122-5p ENSDARG00000000442  slc39a13
	my $key=$tmp[0]."_".$tmp[1];
	if(exists $targetList{$key}){
		my @tmp2=(split/\t/, $targetList{$key}); # dre-miR-122-5p ENSDARG00000000442 slc39a13 hybrid 1
		$tmp2[3]=$tmp2[3].",".$name2;
		$tmp2[4]+=1; # dre-miR-122-5p ENSDARG00000000442 slc39a13 hybrid,mirand 2
		my $value=$tmp2[0]."\t".$tmp2[1]."\t".$tmp2[2]."\t".$tmp2[3]."\t".$tmp2[4];
		$targetList{$key}=$value;
	}else{
		my $value=$_."\t".$name2."\t1";
		$targetList{$key}=$value;
	}
}
close LIST2;


open (LIST3, "$list3") or die "$!\n";
while(<LIST3>){
	chomp;
	my @tmp=(split/\t/, $_); # dre-miR-122-5p ENSDARG00000000442  slc39a13
	my $key=$tmp[0]."_".$tmp[1];
	if(exists $targetList{$key}){
		my @tmp2=(split/\t/, $targetList{$key}); # dre-miR-122-5p ENSDARG00000000442 slc39a13 hybrid 1
		$tmp2[3]=$tmp2[3].",".$name3;
		$tmp2[4]+=1; # dre-miR-122-5p ENSDARG00000000442 slc39a13 hybrid,mirand 2
		my $value=$tmp2[0]."\t".$tmp2[1]."\t".$tmp2[2]."\t".$tmp2[3]."\t".$tmp2[4];
		$targetList{$key}=$value;
	}else{
		my $value=$_."\t".$name3."\t1";
		$targetList{$key}=$value;
	}
}
close LIST3;


open (OUT, "> $output") or die "$!\n";
my @all_keys=sort{$targetList{$a} cmp $targetList{$b}} keys %targetList;
foreach (@all_keys){
	my @tmp=(split/\t/, $targetList{$_}); # dre-miR-122-5p ENSDARG00000000442 slc39a13 hybrid,mirand 2
	if($minoverlap<=$tmp[4]){
		print OUT "$targetList{$_}\n";
	}
}
close OUT;


