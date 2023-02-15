##########################    In script documentation	##########################
=head

Created by Iyanar VETRIVEL in March 2014.
Modified by Iyanar VETRIVEL in April 2014 to change the input and output method.

Description : This script reads chain seperated dssp files and extracts the secondary structure and the solvent accessibility information.Output is written to PDB_chainID.dssp.ss and PDB_chainID.dssp.acc files respectively.

Usage : example- perl dssp_to_ss_acc.pl path/input_file path/output_directory.

=cut
##################################################################################

#!/usr/bin/perl
use File::Basename;

chomp($ARGV[0]);#input dssp file
chomp($ARGV[1]);#output directory where you want the PDB_chainID.dssp.ss,PDB_chainID.dssp.acc file to be created.
($name,$path,$suffix) = fileparse($ARGV[0],qr/\.[^.]*$/);

unless(open(FIN,"<","$ARGV[0]"))
{
	print"Cannot open\"$name$suffix\"\n";
	exit;
}
unless(open(FOUT1,">","$ARGV[1]/$name$suffix.ss"))
{
	print"Cannot open \"$name$suffix.ss\"\n";
	exit;
}
unless(open(FOUT2,">","$ARGV[1]/$name$suffix.acc"))
{
	print"Cannot open \"$name$suffix.acc\"\n";
	exit;
}

@dssp_dump=<FIN>;
foreach $line (@dssp_dump)
{
	unless(($line=~ m/\.$/)||($line=~ m/^(\s*#)/))
	{
		push(@dssp_dump_sans_header,$line);
	}
}

#This for loop replaces a chain break with "X" in the secondary structure and "NA" in accesibility 
for($i=0;$i<@dssp_dump_sans_header;$i++)
{
	@elements_of_line=split("",$dssp_dump_sans_header[$i]);	
	if ($elements_of_line[13]!~ m/!/)
	{
		push(@ss,$elements_of_line[16]);
		$acc_temp=join("",($elements_of_line[35],$elements_of_line[36],$elements_of_line[37],$elements_of_line[38]));
		push(@acc,$acc_temp);
	}
	else
	{
		@elements_of_prev_line=split("",$dssp_dump_sans_header[$i-1]);	
		$prev_res_no=join("",($elements_of_prev_line[7],$elements_of_prev_line[8],$elements_of_prev_line[9]));				
		@elements_of_next_line=split("",$dssp_dump_sans_header[$i+1]);	
		$next_res_no=join("",($elements_of_next_line[7],$elements_of_next_line[8],$elements_of_next_line[9]));
		$res_no_diff="$next_res_no"-"$prev_res_no"-1;
		$x= 'X'x $res_no_diff;				
		$NA= ' NA'x $res_no_diff;				
		push(@ss,"$x");
		push(@acc,"$NA");
		push(@acc," ");
	}
}
=head    #This foreach loop replaces a chain break with "!"
foreach $line2 (@dssp_dump_sans_header)
{

	@elements_of_line2=split("",$line2);	
	if ($elements_of_line2[13]!~ m/!/)
	{
	push(@ss,$elements_of_line2[16]);
	$acc_temp=join("",($elements_of_line2[35],$elements_of_line2[36],$elements_of_line2[37],$elements_of_line2[38]));
	push(@acc,$acc_temp);
	}
	else
	{
	push(@ss,"!");
	$acc_temp=join("",($elements_of_line2[35],$elements_of_line2[36],$elements_of_line2[37],$elements_of_line2[38]));
	$acc_temp=~ s/0/!/;
	push(@acc,$acc_temp);

	}
}
=cut
$SS=join("",@ss);
$SS=~ s/\s/C/g;
print FOUT1"$SS\n";
$ACC=join("",@acc);
print FOUT2"$ACC\n";
