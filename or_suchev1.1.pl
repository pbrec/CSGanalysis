#!/usr/bin/perl -w
use strict;

use File::Copy;
use DBI;

#feste werte
my$wd="/Users/philippbrand/Studium/gradSchool/lab/data/transcriptomics/flaImp/assemblies";
my$orfa="queries/Euglossa_Apis_OR_AA_query.fa";
my$ordb="db/OR/Euglossa_Apis_OR_AA_query";

#blast werte
my$eval=0.000001;
my$threads=4;

my($assem_name,$assem_db,$e,$d,$animal);
my$v=1;

#allowed options

my%opts=(
	-h => 1,
	-a => "",
	-v => 1,
	-e => "",
	-d => "",
	-t => ""
);

&help(), unless(@ARGV>=4);

&get_opt();

#change into wdir and create special assem folder
chdir($wd);
umask(0000);
mkdir($assem_name,0777)|| &help("folder $assem_name still exists: $!\n");
mkdir("/Applications/XAMPP/xamppfiles/temp/$assem_name",0777)||&help("folder /Applications/XAMPP/xamppfiles/temp/$assem_name still exists: $!\n");


#connect to mysql 
my@ary;
my$sth;
my$dbh;
my$db ="CSGanalysis";
my$user = "root";
my$pass = "";
my$host="localhost";
$ENV{MYSQL_UNIX_PORT}= '/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock' ;
$dbh = DBI->connect("DBI:mysql:$db:$host", $user, $pass)
	or die "could not connect to server : $DBI::err ($DBI::errstr)\n";

my$ass=&get_ass($opts{'-a'});

my$asID=&parse_assembly($opts{'-a'},$ass);			#includes assembly into mysql db and returns assembly id
my$ch_ass=&save_ass($asID,1,$assem_name,1);		#header of the assembly contigs are changed and saved to assembly/assembly.fa
#make assembly read blastdb! it's stable
&mkbdb($ch_ass,$assem_db,$asID,1,1);						#blastdb from assembly/assembly.fa is created





my$loop=1;
MAIN:{
	print"Loop $loop started\n";

	my$tbn_out = "$assem_name/OR_tblastn_loop$loop.txt";	
	my$bp_out = "$assem_name/ORF_blastp_loop$loop.txt";
	my$orDB = "db/apidae_ORDB/$loop";

	&blast($orfa,$assem_db,$tbn_out,1,$asID,$loop);				#tblastn of ORs vs assembly database

	#&parse_blast($tbn_out,$asID,$loop,'tbn');	#tblastn is included in mysql db
	&import_sql($asID,"tbn",$tbn_out,$loop,3); #my($asID,$bid,$file,$loop)
	
	my$or_hit_out = &save_ass($asID,$loop,$assem_name,2); #get all sequences with tblastn hits

	my$orf_prot_csv = &orf_finder($or_hit_out,$assem_name,$loop,$asID); #find orfs >=150 AS
	
	die"orf-finder wasn't run successfully\n",unless(-e $orf_prot_csv);

					#include ORFs in mysql DB
	&import_sql($orf_prot_csv,$asID,$loop,4);

	my$orf_protout = &save_ass($asID,$loop,$assem_name,3);		#make ORF-fasta


	#blastp of the ORFs vs the OR DB-> has to be created;
	&mkbdb($orfa,$orDB,$asID,$loop,0);							#make OR blastDB
	&blast($orf_protout,$orDB,$bp_out,0,$asID,$loop);			#blastp ORFs vs OR -> identifies false positives
	#&parse_blast($bp_out,$asID,$loop,'bp');		#save blastp to mysqlDB
	&import_sql($asID,"bp",$bp_out,$loop,3);

	&ident_fp($asID,$loop);

#vll noch ne tabelle mit allen ausgeführten cmds????
	my$next_step=0;

	unless($loop==1){
		#DER GROßE VERGLEICH!!!!!!!

		my$new = &compare_lastloop_with_others($loop,$asID,2);
		$next_step=1,if$new;
	}
	else{
		&compare_lastloop_with_others($loop,$asID,1);
	}


	my$clean_or_set = &save_ass($asID,$loop,$assem_name,4),if($loop==1 || $next_step);

	$orfa=$clean_or_set;

$loop++;

redo MAIN, if($loop==2 || $next_step==1);
}


my$resulting_OR_set = &save_ass($asID,$loop-1,$assem_name,5);


$dbh->disconnect();
print"searching for ORs finished after ", $loop-1, " loops\nresults written to $resulting_OR_set\n";

chmod(0755,$assem_name);
`mv run.log $assem_name`,if(-e "run.log");

exit;

##################################################################################

sub scmd{
	my($asID,$loop,$cmd)=@_;
	$cmd=$dbh->quote($cmd);
	my$query=qq{
		INSERT OR_pipe_cmds (asID,lp,cmd)
            	VALUES ($asID,$loop,$cmd)
	};
	$sth = $dbh->do($query)
                or die "Can't do $query: $dbh->errstr\n";
}

sub compare_lastloop_with_others{

	my($loop,$asID,$mode)=@_;
	my@q;
	my$query;

	if($mode==2){
		#the first query removes everything that is found again in a new loop
		$query=qq{
					UPDATE OR_ORFs as ors,
					(SELECT locus,transfrag,seq FROM OR_ORFs where lp!=$loop and asID=$asID) st
					SET ors.uni=0 
					WHERE asID=$asID AND st.locus=ors.locus AND st.transfrag=ors.transfrag AND st.seq=ors.seq AND lp=$loop

		};
		push@q,$query;
		#the second removes OR ORFs that are already present in the library
		$query=qq{
					UPDATE OR_ORFs as ors,
					(SELECT concat(locus,transfrag) as loc,seq,asID, uni,lp FROM OR_ORFs where uni=1 and asID=$asID and lp!=$loop) st
					SET ors.uni=0 
					WHERE locate(ors.seq,st.seq)>=1 AND st.asID=ors.asID and ors.uni=1 and ors.lp=$loop
		};
		push@q,$query;
		#the third second removes sequences from other loops that fit into this loop
		$query=qq{
					UPDATE OR_ORFs as ors,
					(SELECT concat(locus,transfrag) as loc,seq,asID, uni,lp FROM OR_ORFs where uni=1 and asID=$asID and lp=$loop) st
					SET ors.uni=0 
					WHERE locate(ors.seq,st.seq)>=1 AND st.asID=ors.asID and ors.uni=1 and ors.lp!=$loop;
		};
		push@q,$query;
	}
	elsif($mode==1){
		#this query is used in the first loop do remove doubles from the OR ORF set
		$query=qq{
				UPDATE OR_ORFs as ors,
					(SELECT concat(locus,transfrag) as loc,seq,asID, uni,lp FROM OR_ORFs where uni=1 and asID=$asID and lp=$loop group by seq) st
					SET ors.uni=0 
					WHERE locate(ors.seq,st.seq)>=1 AND st.loc!=concat(ors.locus,ors.transfrag) AND st.asID=ors.asID and ors.uni=1 and ors.lp=$loop

		};
		push@q,$query;
	}
	foreach my$query(@q){
		$sth = $dbh->do($query)
        	        or die "Can't do $query: $dbh->errstr\n";
    }

    my$new=&retrieve_sql($asID,$loop,2);
    &scmd($asID,$loop,$query);
    return $new;
}

sub ident_fp{

	my($asID,$loop)=@_;

	my$query=qq{
		UPDATE OR_ORFs as ors,
		(SELECT locus as bpC, transfrag as bpT, orfnum as bpnum
		FROM OR_blastp where asID=$asID and lp=$loop and evalue<=1e-6 
		group by locus, transfrag,orfnum) st
		SET ors.tp=1 
		WHERE st.bpC=ors.locus AND st.bpT=ors.transfrag AND st.bpnum=ors.orfnum AND asID=$asID AND lp=$loop
	};
	$sth = $dbh->do($query)
                or die "Can't do $query: $dbh->errstr\n";
    &scmd($asID,$loop,$query);

}


sub orf_finder{

	my($in,$ass_name,$loop,$asID)=@_;
	my$DNAout="$ass_name/orf_DNAout_loop$loop.fa";
	my$Protout="$ass_name/orf_Protout_loop$loop.tsv";

	my$cmd="./orf-finder.pl $in $DNAout $Protout";
	print"CMD: $cmd\n";
	system($cmd);
	print"DONE\n";

	&scmd($asID,$loop,$cmd);

	return$Protout;

}

sub save_ass{

	my$mode=pop@_;
	my$outfile=pop@_;
	my($asID,$loop)=@_;
	my$query;
	if($mode==1){
		#changes assembly header outfile to >locus_transfrag

		
		$outfile="$outfile/$outfile.fa";
		print "change header\n";
		$query=qq{select concat(">",locus,"_",transfrag,"\n",seq) into outfile '/Applications/XAMPP/xamppfiles/temp/$outfile' fields terminated by "" escaped by "" from assemblies where asID=$asID} ;
		print "Done\n";

	}
	elsif($mode==2){

		#creates fasta-file with OR-tbn hits <=1e-6
		
		
		$outfile.="/tbn_or_hits_loop$loop.fa";
		#soap-evir-1-1-merged/tbn_or_hits_1.fa
		$query=qq{select concat(">",st.locus,"_",st.transfrag,"\n",ass.seq) 
					into outfile '/Applications/XAMPP/xamppfiles/temp/$outfile' 
					fields terminated by "" escaped by "" 
					from 
					assemblies as ass, 
						(select locus,transfrag 
							from 
							OR_tblastn 
							where 
							evalue<=1e-6 and asID=$asID and lp=$loop group by locus, transfrag) st 
					where 
					asID=$asID and st.locus=ass.locus and st.transfrag=ass.transfrag};
	}
	elsif($mode==3){ #saves AS-Seqs with OR ORFs in fa format HEADER >loc_tr_orf_num
		

		$outfile.="/ORF_loop$loop.prot.fa";
		$query=qq{
			select concat(">",locus,"_",transfrag,"_",orfnum,"\n",seq) 
			into outfile '/Applications/XAMPP/xamppfiles/temp/$outfile' 
			fields terminated by "" escaped by "" from OR_ORFs where asID=$asID and lp=$loop

				} ;

	}
	elsif($mode==4){
		
		$outfile.="/clean_ORF-set_loop$loop.prot.fa";
		$query=qq{
				SELECT concat(">",orfnum,"\n",seq)
				into outfile '/Applications/XAMPP/xamppfiles/temp/$outfile'
				fields terminated by "" escaped by "" 
				from OR_ORFs 
				where asID=$asID and lp=$loop and tp=1 and uni=1
				};
	}
	elsif($mode==5){
	
		$outfile.="/resulting_OR_set.fa";
		$query=qq{
				SELECT concat(">",locus,"_",transfrag,"_",orfnum,"\n",seq)
				into outfile '/Applications/XAMPP/xamppfiles/temp/$outfile'
				fields terminated by "" escaped by "" 
				from OR_ORFs 
				where asID=$asID and tp=1 and uni=1
				};
	}
	$sth = $dbh->do($query)
                or die "Can't do $query: $dbh->errstr\n";
    		
    print"DONE\n";
	move("/Applications/XAMPP/xamppfiles/temp/$outfile", "/Users/philippbrand/Studium/gradSchool/lab/data/transcriptomics/flaImp/assemblies/$outfile")|| die "cannot move file: $!\n";
    &scmd($asID,$loop,$query);


    return($outfile);

}

sub parse_blast{
	#needs blastfile assembly-id loop-nr blasttype
	print"import blast\n";
	my($open,$asID,$loop,$type)=@_;
	

	open(IN,$open)|| die "Cannot open $open: $!\n";

	while (my$line=<IN>){
		chomp($line);
		my@split=split /\s+/, $line;
		unshift(@split,$type);
		unshift(@split,$asID);
		push(@split,$loop);

		&import_sql(@split,3);
	}

	#print"DONE\n";
}

sub parse_assembly{
        #the assembly fa file is being parsed and introduced to the OR_blast database
        print"Import assembly into database\n";
        my($path,$ass)=@_;
        #my$ass=&get_ass($path);
        open(IN,$path)|| die "cannot open $path to parse it:$!\n";
        my%pa;
        my$head=0;
        while (my$line=<IN>){
                chomp($line);
                if($line=~/^>/){

                        if($ass eq "s"){
                                my@parse=split /\s+/,$line;
                                #$parse[0]=~s/>C//;
                                $head+=1;
                    			$parse[0]=$head;

                                splice(@parse,1,0,1);#sql format: Locus transcr extra
                    
              
                                $pa{$head}=\@parse;
                        }
                        elsif($ass eq "t"){
                        	$head+=1;
                        	my@parse=split /\s+/,$line;
                        	my$locraw=shift@parse;
                        	$locraw=~/c(\d+)_g(\d+)_i(\d+)/;
                        	my$locus="$1" . "0" . "$2";
                        	my$trans="0" . "$3";
                        	my$extra=join " ",@parse;
                        	@parse=("$locus","$trans","$extra");

                        	$pa{$head}=\@parse;
                        }
                        elsif($ass eq "o"){
                        	$head+=1;
                        	$line=~/Locus_(\d+)_Transcript_(\d+)\/\d+_Confidence_(\d+\.?\d*)/;
                        	my@save=("$1","$2","$3");
                        	$pa{$head}=\@save;
                        }
                }

                else{
                        $line=~s/\s+//g;
                        ${$pa{$head}}[3].=$line;
                }

        }
        close IN;
        &import_sql($opts{'-t'},$ass,$opts{'-e'},$opts{'-d'},1);

        my$asID=&retrieve_sql($opts{'-t'},$ass,$opts{'-e'},$opts{'-d'},1);

		foreach my$data(keys %pa){
        	
        	&import_sql($asID,$opts{'-t'},${$pa{$data}}[0],${$pa{$data}}[1],${$pa{$data}}[2],length(${$pa{$data}}[3]),${$pa{$data}}[3],2);
    	}

    	&basic_assem_stats($path,$asID);


    	print"DONE\n";
    	return $asID;
}

sub basic_assem_stats{
	my($assembly,$asID)=@_;
	my$cmd="./basic_statistics.pl $assembly";
	my$results=`$cmd`;
	die "stat_retrieval did not work\n",unless($results or $results=~/^\d+/);
	my@results=split /\s+/,$results;
	unshift@results,$asID; 
	push@results,5;

	&import_sql(@results);
	&scmd($asID,1,$cmd);
}

sub retrieve_sql{
	my$mode=pop@_;
	my$query;
	if($mode==1){
		my($anim,$ass,$e,$d)=@_;

		$query=qq/SELECT id FROM assembler_set WHERE animal="$anim" AND assembler="$ass" AND e=$e AND d=$d/;
		
	}
	elsif($mode==2){
		my($asID,$loop)=@_;
		$query=qq/
				SELECT count(locus) 
				FROM OR_ORFs 
				WHERE tp=1 AND uni=1 AND asID=$asID AND lp=$loop
		/;
	}
	else{&help("Retrieval mode for sql unknown\n")}

	$sth  = $dbh->prepare($query)
                or die "Can't prepare $query: $dbh->errstr\n";

    $sth->execute
                or die "can't execute the query: $sth->errstr";
    my$val;

    while(@ary=$sth->fetchrow_array()){
    	#if($mode==1){
    		$val=shift@ary;
    	#}
    }
	&fin();
    return($val);
}

sub import_sql{
	my$mode=pop(@_);
	my$query;
	my$prefix="/Users/philippbrand/Studium/gradSchool/lab/data/transcriptomics/flaImp/assemblies/";
	if($mode==1){
		my($anim,$ass,$e,$d)=@_;
		$query=qq/
            	INSERT assembler_set (animal,assembler,e,d)
            	VALUES ("$anim","$ass",$e,$d)
               /;
        
    }
    elsif($mode==2){
    	my($id,$animal,$loc,$tf,$extra,$lens,$seq)=@_;
    	$extra=$dbh->quote($extra);
    	$query=qq/
                INSERT assemblies (asID,animal,locus,transfrag,extra,lens,seq) 
                VALUES ($id,"$animal",$loc,$tf,$extra,$lens,"$seq")
               /;
    }
    elsif($mode==3){ #imports blasts in mysql
    	#my($asID,$bid,$qid,$sid,$lens,$qs,$qe,$ss,$se,$eval,$bit,$pid,$nid,$loop)=@_;

    	my($asID,$bid,$file,$loop)=@_;
    	my$orf_num=0;		#wenn 0 gibt es keine!
    	my$transfrag=undef;
    	

    	if ($bid eq "tbn"){

    			$query=qq/
                	LOAD DATA LOCAL INFILE '$prefix$file' IGNORE INTO TABLE OR_tblastn FIELDS TERMINATED BY '\t' (qseqid,\@hit,lens,qstart,qend,sstart,send,evalue,bitscore,pident,nident) 
	SET locus=SUBSTRING_INDEX(\@hit,'_',1), transfrag=SUBSTRING_INDEX(\@hit,'_',-1), asID = $asID, lp=$loop
               /;

    	}
    	elsif($bid eq "bp"){
    		

    		$query=qq/
    		LOAD DATA LOCAL INFILE '$prefix$file' IGNORE INTO TABLE OR_blastp FIELDS TERMINATED BY '\t' (\@hit,sseqid,lens,qstart,qend,sstart,send,evalue,bitscore,pident,nident) 
			SET transfrag=substring_INDEX(substring_INDEX(\@hit,'_',2),'_',-1), locus=SUBSTRING_INDEX(\@hit,'_',1), orfnum=SUBSTRING_INDEX(\@hit,'_',-1), asID = $asID, lp=$loop
               /;
    	}

    }
    elsif($mode==4){ #imports ORFs into mysql DB
    	my($file,$asID,$loop)=@_;
    	$query=qq/
    			LOAD DATA LOCAL INFILE '$prefix$file' IGNORE INTO TABLE OR_ORFs FIELDS TERMINATED BY '\t' (locus,transfrag,orfnum,\@seq,dnaseq) 
				SET seq=\@seq, lens=length(\@seq), asID = $asID, lp=$loop
    		/;
    }
    elsif($mode==5){
    	my($asID,$N25,$N50,$N75,$N95,$N10,$N20,$N30,$N40,$N60,$N70,$N80,$N90,$longest,$shortest,$mean,$median,$bp)=@_;
   		$query=qq/
    			UPDATE assembler_set 
    			SET N25=$N25,N50=$N50,N75=$N75,N95=$N95,N10=$N10,N20=$N20,N30=$N30,N40=$N40,N60=$N60,N70=$N70,N80=$N80,N90=$N90,longest=$longest,shortest=$shortest,mean=$mean,median=$median,bp=$bp
    			WHERE id=$asID
    			/;
    }
    else{&help("Could not locate sql import mode\n");}

    $sth = $dbh->do($query)
                or die "Can't do $query: $dbh->errstr\n";
    

}

sub fin{
	$sth->finish();
}

sub get_ass{
	#finds the assembler by looking at the fasta head
	my$path=shift;
	my$ass;
	my@head=`head $path`;
        foreach my$line(@head){
                next,unless($line=~/^>/);
                if($line=~/^>C\d+/){
                        $ass="s";
                }
                elsif($line=~/^>Locus_\d+_Transcript/){
                        $ass="o";
                }
                elsif($line=~/^>c\d+_/){
                        $ass="t";
                }
        }
        die"don't know assembly format\n", unless($ass);
        return($ass);
}

sub mkbdb{

        my($in,$out,$asID,$loop,$dbtype)=@_;
        $dbtype=$dbtype ? "nucl" : "prot";

        my$makedb="makeblastdb -in $in -out $out -dbtype '$dbtype'";
        print "CMD: $makedb", if($v);
        system($makedb);
        print "DONE\n",if($v);
        &scmd($asID,$loop,$makedb);
}

sub get_opt{
	my$c=-1;
	foreach my $opt (@ARGV){
		$c++;
        next,unless($opt=~/^-/);
        &help(), if($opt eq '-h');

        &help("param missing behind $opt\n"), if($ARGV[$c+1]=~/^-/);

        if (exists $opts{"$opt"}){
        	$opts{$opt}=$ARGV[$c+1];
       		if($opt=~/-e|-d|-v/){
       			&help("argument for opt $opt must be an int\n"),unless($ARGV[$c+1]=~/^\d+\.?\d*$/);
       		} 	
       		elsif($opt=~/-t/){
       			&help("argument for opt $opt must be v or d\n"),unless($ARGV[$c+1]=~/^v$|^d$|^b$|^i$|^f$/);
       		}
       		elsif($opt=~/-a/){
       			&help("file $ARGV[$c+1] does not exist\n"),unless(-e $ARGV[$c+1]);
       			$assem_name=$opts{$opt};
                $assem_name=~s/^[\w\W]+\///;
                #$assem_name=~s/\.[\w\W]+$//;
                $assem_name=~s/\.fasta$//;
                $assem_db="db/$assem_name/$assem_name";
       		}
        }
        else{
        	&help("unknown option: $opt\n");
        }
        
	}

}
sub blast{
        my($q,$db,$o,$blast,$asID,$loop)=@_;

        $blast= $blast? "tblastn" : "blastp";

        my$blastex="$blast -query $q -db $db -out $o -outfmt '6 qseqid sseqid length qstart qend sstart send evalue bitscore pident nident' -evalue $eval -num_threads $threads";
        print"CMD:$blastex\n";
        system($blastex);
        print"DONE\n";
        &scmd($asID,$loop,$blastex);

}

sub help{
	my$warn=shift;
	$warn=$warn?$warn:" ";
	my$usage=<<_USAGE_;

options:
	
	required:
        
        -a assembly.fa
        -e cov_cutoff
        -d 2nd setting
        -t animal: v d

    misc:

        -h help
        -v verbose: 1/0 default 1
_USAGE_
;
	
	die $warn , $usage;

}

