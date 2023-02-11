##########################    In script documentation	##########################
=head

Created by Iyanar VETRIVEL in March 2014.
Modified by Iyanar VETRIVEL in April 2014 to change the input and output methods.

Description : This script reads chain seperated pdb files and extracts the amino acid sequence based on the ATOM records. Output is written to PDB_chainID.pdb.aa.

Usage : example- perl pdb_to_aa.pl path/input_file path/output_directory.

=cut
##################################################################################

#!/usr/bin/perl

use File::Basename;
chomp($ARGV[0]);#input pdb file
chomp($ARGV[1]);#output directory where you want the PDB_chainID.pdb.aa file to be created.

##############################

$flag=0;

###############################




($name,$path,$suffix)= fileparse($ARGV[0],qr/\.[^.]*$/);

unless(open(FIN1,"$ARGV[0]"))
{
	print"Cannot open \"$name$suffix\"\n";
	exit;
}

unless(open(FOUT,">","$ARGV[1]/$name$suffix.aa"))
{
	print"Cannot create \"$name$suffix.aa\"\n";
	exit;
}



@elements_of_filename1=split(/\./,$name);
@elements_of_filename2=split(/\_/,$elements_of_filename1[0]);
$chain=$elements_of_filename2[1];




@pdb_dump=<FIN1>;
foreach(@pdb_dump)
{
	if($_=~ /^TER/)
	{
		$last_res_num=substr($_,22,4);
	}
}
$last_res_num=~ s/\s+//g;

foreach $p(@pdb_dump)
{
if($p=~ /^TER/)
{
	last;
}



##################################################

elsif($p=~ /^REMARK 470.{9}$chain\s+$last_res_num[(.+\bO\b.+)|(.+\bCA\b.+)|(.+\bC\b.+)|(.+\bN\b.+)]/)
{
#print"$p\n";
$flag=1;
}



##################################################

elsif(($p=~/^(\bATOM\b)|^(\bHETATM\b)/) && (substr($p,16,1)=~/[A| ]/))
		{
			push(@atom_records,$p);
		}	
}

foreach $k(@atom_records)
{
	
	$three_letter_code= substr($k,17,3);
	$three_letter_code=~ s/\s+//g;
	$res_number= substr($k,22,5);
	$res_number=~ s/\s+//g;
	unless($three_letter_code eq 'HOH')
		{
			push(@non_unique_aa_res,"$three_letter_code,$res_number");		
		}
	
}

%seen1=();
@unique_aa_res = grep { ! $seen1{$_}++ } @non_unique_aa_res;

#If $flag is equal to 1 then it means there are missing atoms in the last residue, hence exclude it from the array @unique_aa_res by using pop.
if($flag=='1')
{pop(@unique_aa_res);}


#This section of the code generates an array @missing_res, which contains the "missing residue,residue number" from the REMARK 465 section of the PDB header 
#####################################################################################


for($i=0;$i<$#pdb_dump;$i++)
{
	if($pdb_dump[$i]=~/REMARK 465   M RES C SSSEQI/)
	{
		for($k=($i+1);$k<$#pdb_dump;$k++)
		{
			if($pdb_dump[$k]=~/^(\bREMARK 465\b)/)
			{
				@elements_of_missing=split(/\s+/,$pdb_dump[$k]);
				chomp($chain);				
				chomp($elements_of_missing[3]);
				if($elements_of_missing[3] eq $chain)
				{
					unless($elements_of_missing[4]>$last_res_num)
					{
					$res_comma_res_no=join(",",($elements_of_missing[2],$elements_of_missing[4]));
					$res_comma_res_no=~ s/\s+//g;					
					push(@missing_res,$res_comma_res_no);
					}

				}

			}

		}
	}
}

####################################################################################



#print"\n\nMissing residues-@missing_res\n\n";
#print"Chain-$chain\n\nunique_aa_res before appending\n@unique_aa_res\n";
for($l=0;$l<$#unique_aa_res;$l++)
{
	@elements_of_unique_aa_res=split(/,/,$unique_aa_res[$l]);
	@next_elements_of_unique_aa_res=split(/,/,$unique_aa_res[$l+1]);
	chomp($elements_of_unique_aa_res[1]);		
	if($elements_of_unique_aa_res[1]=~m/[A-Z]$/)
	{
		@parts_of_elements_of_unique_aa_res_1=split("",$elements_of_unique_aa_res[1]);
		$m=pop(@parts_of_elements_of_unique_aa_res_1);
		++$m;
		push(@parts_of_elements_of_unique_aa_res_1,$m);		
		$plus_one=join("",@parts_of_elements_of_unique_aa_res_1);
		$plus_one=~ s/\s+//g;		
		#print"$plus_one\n\n\n";
	}	
	else
	{$plus_one="$elements_of_unique_aa_res[1]"+1;}
	if($next_elements_of_unique_aa_res[1]ne$plus_one)
	{
		unless($next_elements_of_unique_aa_res[1]=~ m/$elements_of_unique_aa_res[1]A/)
		{
		$counter=@missing_res;
		#print"$counter\n";
		for($n=0;$n<$counter;$n++)
		{
			@elements_of_missing_res=split(/,/,$missing_res[$n]);
			if(($elements_of_missing_res[1] eq $plus_one)||($elements_of_missing_res[1] eq "$elements_of_unique_aa_res[1]"+1))
			{
				splice(@unique_aa_res,($l+1),0,$missing_res[$n]);
				$n++;				
			}
		}
		}
	}

}
#print"Chain-$chain\n\nunique_aa_res after appending\n@unique_aa_res\n";
foreach(@unique_aa_res)
	{
		@new_elements_of_unique_aa_res=split(/,/,$_);
		push(@three_letter_aaseq,$new_elements_of_unique_aa_res[0]);

	}
#print"\n\nThree letter aaseq:-\n@three_letter_aaseq\n";

%aa=("ALA"=>"A","ASX"=>"B","CYS"=>"C","ASP"=>"D","GLU"=>"E","PHE"=>"F","GLY"=>"G","HIS"=>"H","ILE"=>"I","LYS"=>"K","LEU"=>"L",
"MET"=>"M","MSE"=>"M","ASN"=>"N","PRO"=>"P","GLN"=>"Q","ARG"=>"R","SER"=>"S","THR"=>"T","VAL"=>"V","UNK"=>"U","TRP"=>"W","TYR"=>"Y","GLX"=>"Z");

for($p=0;$p<=$#three_letter_aaseq;$p++)
{

	if(exists($aa{$three_letter_aaseq[$p]}))
	{
		$single_letter_aaseq[$p]=$aa{$three_letter_aaseq[$p]};
	}
	else
	{
		$single_letter_aaseq[$p]="X";
	}
	
}
$joined_single_letter_aaseq=join("",@single_letter_aaseq);
print FOUT "$joined_single_letter_aaseq";
