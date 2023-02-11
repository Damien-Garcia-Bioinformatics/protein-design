##########################    In script documentation	##########################
=head

Created by Iyanar VETRIVEL in February 2014.
Modified by Iyanar VETRIVEL in April 2014 to change the input and output methods.

Description : This script reads DSSP files and separates them into their individual chains. Output is written as PDB_chainID1.dssp, PDB_chainID2.dssp, PDB_chainID3.dssp..... 		      

Usage : example- perl dssp_separator.pl path/input_file path/output_directory.

=cut
##################################################################################

#!/usr/bin/perl

use File::Basename;
chomp($ARGV[0]);#input dssp file
chomp($ARGV[1]);#output directory where you want the PDB_chainID.dssp file to be created.
($name,$path,$suffix)= fileparse($ARGV[0],qr/\.[^.]*$/);


unless(open(FIN,"<","$ARGV[0]"))
{
	print"Cannot open \"$name$suffix\"\n";
	exit;
}

print"Splitting the dssp file....\n";
@dssp_dump=<FIN>;
$full_variable=join("",@dssp_dump);
@parts=split(/(\n.+#.+\n)/,$full_variable);
$header=join("",$parts[0],$parts[1]);
@dssp_lines=split(/\n.+!\*.+\n/,$parts[2]);
for($i=0;$i<$#dssp_lines+1;$i++)
{
	@dssp_parts=split("",$dssp_lines[$i]);	
	$chain_id=$dssp_parts[11];
	open(FOUT,">","$ARGV[1]/$name\_$chain_id.dssp");
	$full_data=join("",$header,$dssp_lines[$i]);
	print FOUT "$full_data";
}

