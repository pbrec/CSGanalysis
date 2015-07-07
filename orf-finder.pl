#!/usr/bin/perl -w
use strict;

my$lens=150; #sets the minimum length of ORFs to be reported

die "This version needs 3 arguments!\n 1. Input: (batch) DNA sequences in fastaformat\n 2. DNAout file (is created by the script)\n 3. Protout file (also created)\n",unless($ARGV[2]); 

my$open=$ARGV[0];	#datei mit DNAsequenzen
open(IN,$open)||die"sorry, cannot open $open: $!\n";

my%fa;
my$hold;
my$locus=0;
#speichert fasta file nach oases loci geordnet
print"start to open the fasta file\n";
while(my$line=<IN>){
	chomp($line);

	if($line=~/^>/){
			my@locus=split /_/,$line;
			$locus=$locus[0];
		${$fa{$locus}}{$line}=();
		$hold=$line;
	}
	else{
		${$fa{$locus}}{$hold}.=$line;
	}
}
close IN;
print "opening completed\n";


#orf finder
print "orf-finder has started!\n";
my%orf;


foreach my$loc(keys %fa){
	#für jeden locus
	my$hits=0;
	my%tmp;
	foreach my$ctg(keys %{$fa{$loc}}){
		#in jedem Ctg
	
		my$seqa=${$fa{$loc}}{$ctg};
		#find orfs mit regex
		#print"ctg: $ctg\n";

		my$dir="fwd";
		my($fwd,$rev)=(0,0);


		for(my$i=0;$i<2;$i++){		#2 mal: fwd und rev

			my$seq_b=$seqa;
			if($i){
				#make rc
				$seq_b=reverse($seq_b);
				$seq_b=~tr/AGTC/TCAG/;
				$dir="rev";
			}

			
			#print "direction: $dir\n";

			for(my$j=0;$j<3;$j++){
				#for all 3 frames
				$seq_b=~s/\w{1,1}//, if($j);
				
				my$seq=$seq_b;
				my$frame=$j+1;
				#print"Frame: $frame \n";
				while( $seq=~/(?:...)*?(?<hit>ATG(?:...)*?(?:TAA|TGA|TAG))/i ){
					
					if(length($+{hit}) >= $lens){
						
						$hits++;
						#print "hit$hits\nlens:" . length($+{hit}) . "\n$+{hit}\n";

						${$tmp{$+{hit}}}{$ctg}=();
					}
					$seq=$';
				}
			}
			#$dir eq "fwd" ? $fwd=$hits : $rev=$hits-$fwd;
		}
		
		

	}
	#print"insgesamt $hits orfs gefunden\n";

	&reduction(\%tmp);

	#print"Non reduced ORFs\n";

	#my$a=1;
	#foreach my$hits(keys %tmp){

		
	#	print ">$a\n$hits\n";
	#	$a++;
	#}
}

#exit(0);
#print"Die Hits:\n";
open(OUT,">$ARGV[1]");
foreach my$hits(keys %orf){

	foreach my$names(keys %{$orf{$hits}}){


		print OUT "$names\n";
	}
	print OUT "$hits\n";
}
close OUT;


#TRANSLATE ORFs and print them
#my%translations;

my@date=localtime(time);
my$y=$date[5]+1900;
$date[4]++;
my$date="$date[3]_$date[4]_$y";

#open(OUT1,">uniqueORFs_$date.fasta");
open(OUT1,">$ARGV[2]");

print "translate .";


my$orf_numb=1;

foreach my$unique (keys%orf){

	my@names=keys(%{$orf{$unique}});



	&translatE($unique,$names[0]);

}
close(OUT1);
print"\neverything done, Goodbye\n";


exit(0);
#WENN DAS SKRIPT SCHEISSE LAUFT $add enTfernen!!!!
sub translatE{

	my($dna,$name)=@_;
	my@split=split /_/,$name;
	$name=$split[0];
	$name=~s/^>//;
	my$tr=$split[1];
	#my$add=$split[2]=~/\d+/?$split[2]:"";
	my%codons=("TTT" => "F","TCT" => "S","TAT" => "Y","TGT" => "C","TTC" => "F","TCC" => "S","TAC" => "Y","TGC" => "C","TTA" => "L","TCA" => "S","TAA" => "", "TGA" => "", "TTG" => "L","TCG" => "S","TAG" => "", "TGG" => "W","CTT" => "L","CCT" => "P","CAT" => "H","CGT" => "R","CTC" => "L","CCC" => "P","CAC" => "H","CGC" => "R","CTA" => "L","CCA" => "P","CAA" => "Q","CGA" => "R","CTG" => "L","CCG" => "P","CAG" => "Q","CGG" => "R","ATT" => "I","ACT" => "T","AAT" => "N","AGT" => "S","ATC" => "I","ACC" => "T","AAC" => "N","AGC" => "S","ATA" => "I","ACA" => "T","AAA" => "K","AGA" => "R","ATG" => "M","ACG" => "T","AAG" => "K","AGG" => "R","GTT" => "V","GCT" => "A","GAT" => "D","GGT" => "G","GTC" => "V","GCC" => "A","GAC" => "D","GGC" => "G","GTA" => "V","GCA" => "A","GAA" => "E","GGA" => "G","GTG" => "V","GCG" => "A","GAG" => "E","GGG" => "G");

	#print "NT:\n$dna\n";
	my@triple=split//,$dna;
	my$as;
	
		while(my@tpl=splice @triple,0,3){
		my$tpl=join("",@tpl);
		$tpl=uc($tpl);
		if ($tpl=~m/N/){
			$as.="X";
		}
		else{	
			unless (exists $codons{$tpl}){
				print "Unknown template: $tpl -> I'll add an X\n";
				$as.="X";
			}
			else{
				$as.=$codons{$tpl};
			}
		}
		#print"\rtranstlate $codons{$tpl}";
	}
	
	print OUT1 "$name\t$tr\t$orf_numb\t$as\t$dna\n";
	$orf_numb++;
}



sub reduction{
	#print"reduction of ORF-Set started\n";
	my%hit=%{$_[0]};
	#die geffundenen Sequenzen werden der länge nach sortiert (kürzeste zuerst)
	my@all_keys=sort {length($a)<=>length($b)}(keys %hit);


REDO:	foreach my$seq (sort {length($a)<=>length($b)}(keys %hit)){
		#print "Length:\n" . length($seq) . "\n";
		shift @all_keys;

		if(@all_keys){
			foreach my$others(@all_keys){
				#if (index($seq,$others) == -1 ){ #geht so nicht! auf fraeshift muss geachtet werden!
				my$test_seq=$seq;
				$test_seq=~s/(?:...)$//;	#das stop codon entfernen!
				
				#print"Seq:\n$seq\nTest_Seq:\n$test_seq\n";
				#exit(0);

				#schauen, ob ein orf teilstring eines anderen ist. Aber natürlich ohne stop codon
				if($others=~/(...)*?$test_seq/){

					delete ($hit{$seq});
					next REDO;
				}
			}
		}

	}

	#print"ORFs after reduction: " . keys(%hit) . "\n";

	foreach my$add(keys %hit){
		$orf{$add}=$hit{$add};
	}
	#print"reduction completed\n";
}
