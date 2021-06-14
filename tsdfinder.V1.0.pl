#! /usr/bin/perl -w

{ # Define the main body of code as a block so that variables
  # defined in main are not global.

# The bl2seq (blast 2 sequences) code must be executable from
# the command line from the machine on which you are running
# this code!  The bl2seq code is available from the NCBI ftp
# site (ftp://ftp.ncbi.nih.gov/blast/). 

# The RepeatMasker program must be run first on the sequence
# you want to analyze for L1s.  The only RepeatMasker output files 
# that are needed for TSDfinder is/are the *.out file(s).  
# ALSO, when running RepeatMasker, use a custom library that
# contains a single L1 reference sequence.  We use the 
# L1.3 sequence.

# This program is written to process chromosomes, one at a time.
# So it assumes that each sequence file contains the 
# complete sequence of a chromosome and thus the corresponding
# RepeatMasker file contains the L1 ananotation for that chromosome.

# New directories will be created in the directory in which
# this program is run.  These directory names will be named 
# according to the chromosome being analyzed.

# TSDfinder must also have access to the original sequence
# files in fasta format that are being analyzed. 
# TSDfinder makes the following assumptions about the sequence:
# 1) The sequence files are zipped (e.g. have the 
#    .gz suffix on the Unix OS) and end with the suffix ".fa".
# 2) The sequence files that correspond
#    to the RepeatMasker *.out files that are being analyzed have
#    the same name (except, of course, for the lack of the ".out"
#    suffix and additional ".gz" suffix on the sequence files).
# 3) The definition line(s) start(s) with ">gi".

# Make sure the *.out file(s) from RepeatMasker
# is/are in the same directory in which you're running TSDfinder.

# MAKE THESE CHANGES BEFORE RUNNING TSDFINDER:
# List the names of the RepeatMasker files on line 140.
# Change line 77 so that it reflects the directories
# in which the sequence being analyzed is located.

# How program works:
# For each *.out file in current directory, extract
# gi and coordinates of "full insertions" (correspond
# to single complete insertion event, as judged by
# coordinates in query and L1.3, as reported by 
# RepeatMasker).  Merge coordinates of segments that 
# constitute the same insertion event (inverted-deleted).
# Mark these with INV flag in MERGED output.

# TSDfinder assigns each L1 element a unique ID. 
# The output files include MERGED which contains the coordinates
# of all the L1s in both the input seuquence and L1 reference
# sequence -- very similar to the original RepeatMasker *.out file.
# INVERSIONS contains the coordinates of all inverted L1s and indicates
# where the breakpoint between the 2 segments is.
# SUMMARY contains the general summary and statistics of TSDfinder results.
# PREINSERTION contains the pre-insertion locus of each L1 for which TSDs
# were found.  It also gives the TSD sequence, the TSD score, the polyA tail 
# (if no tail found but a TSD pair was found, there'll be an 'X'), and it
# indicates whether the L1 is a 3' transduction event.

# This same dir will  contain all of the Chromosome dirs
# with the TSD finding results.

$| = 1;

use strict;

# Allow $insertion_event_id and all constants to be global
use vars qw($insertion_event_id $AllInsertions %AllVsIntact $DELTAGI $DELTAL1 $THR_PR_END $FIVE_PR_FLANK $THREE_PR_FLANK $L1_TIPS $FLANKS_PREINSERTION $PREFIXtoID $DirOfGenomeSeqs);

$DirOfGenomeSeqs = " ";

# Program set-up
my(@Files, $NumFiles, $index, $MkDir, $MkMainDir);
# Variables used in collecting information about the 
# L1s from the *.out files.
my($n, $filename,$chromosome); 
my($line_number, $Line, $line, $new_insertion_events);
my(@gi, @L1left, $count, @PercId, @start, @end);
my(@strand, @L1end, @L1start, @ExtraSeq, @Piece, $PercIden, $Length1, $Length2);
# Checking quality of each L1 reported by RepeatMasker
my($inversionFLAG, $deltagi, $deltaL1);
# Setting up TSD finding
my($ins_event, $no_TSDs_counter, $standard_insertion_counter, $TSDs_wo_ATail);
my($transduction_events_counter, $KeepIt, $ConvertId);
# Arrays of information for L1 insertions that will be analyzed for TSDs
my(@Gi, @Strand, @Start_in_Gi, @End_in_Gi, @L1left_after_merge, @L1_length, @inversionFLAG);

# Files for bl2seq analysis
my($a_file, $b_file, $Out_file, $i_file, $j_file, $o_file);
# Evaluating bl2seq results
my($hspFLAG, $HSPcounter, $TSDscore, $GotAsForEvent, $PolyAScore, $CurrentBestTSD);
my($HSPidentities, $HSPlength, $transduction_state, $transduction_filename);

my($TSD_5pr_start_inGi, $TSD_5pr_start, $Start_in_Gi, $preinsertion_locus_start_inGi);
my($TSD_3pr_end_inGi, $TSD_3pr_end, $End_in_Gi, $preinsertion_locus_end_inGi); 

my($pos_before_5prTSD, $pos_after_3prTSD);
my($transd_segment_start, $transd_segment_end);

my($preins_5pr_sequence, $preins_3pr_sequence);
my($HSP_5pr_start, $HSP_5pr_sequence, $HSP_5pr_end, $TSDstart, $TSDend);
my($HSP_3pr_start, $HSP_3pr_sequence, $HSP_3pr_end, $ThrPrTSDstart, $ThrPrTSDend);
my($TSD_5pr_sequence, $TSD_3pr_sequence, $TSD, @Tsd_5pr, @Tsd_3pr, $tsd_length);
my($bp, $pos_before_3prTSD, $pos_after_5prTSD, @FIVE_PR, @THREE_PR);

# PolyA - related variables
my(@ReturnVals, $Result, $Score, $Motif, $TSDTail, $LengthOfTail);
my($TSDPolyATail, $TSDPolyATailMotif);

# Files for keeping L1 information for each chromosome
my($MergeFile, $InversionFile, $PreInsFile, $StatFile);

# Set parameters for merging L1 segments (for rearranged L1s):
$DELTAGI = 30;	
# Length of indel allowed in L1:
$DELTAL1 = 50;	
# Minimum number of bp required at 3' of L1 to be analyzed for TSDs.
$THR_PR_END = 35;   

# Set parameters for flanking sequence retrieval
$FIVE_PR_FLANK = 100;       # bp	|
$THREE_PR_FLANK = 3000;     # bp	| for TSD finding step
$L1_TIPS = 15;		    # bp 	|
$FLANKS_PREINSERTION = 50;  # bp--reconstruction of the pre-insertion locus

# In this sub, confirm that user wants to use the above default parameters.
&IntroductoryRemarks;

# In the following @Files array, place the names of the RepeatMasker *.out files that
# you want to analyze.  The following is an example:
## @Files = qw(hs_chr1.fa.out hs_chr2.fa.out hs_chr3.fa.out hs_chr4.fa.out hs_chr5.fa.out hs_chr6.fa.out hs_chr7.fa.out hs_chr8.fa.out hs_chr9.fa.out hs_chr10.fa.out hs_chr11.fa.out hs_chr12.fa.out hs_chr13.fa.out hs_chr14.fa.out hs_chr15.fa.out hs_chr16.fa.out hs_chr17.fa.out hs_chr18.fa.out hs_chr19.fa.out hs_chr20.fa.out hs_chr21.fa.out hs_chr22.fa.out hs_chrY.fa.out hs_chrX.fa.out hs_chrUn.fa.out);

@Files = qw();

$NumFiles = @Files;

$n = 0;  # Counter for *.out files in the current directory

# Open each *.out file in current directory; for each *.out file
# create *.out.merged output

for ($index = 0; $index < $NumFiles; $index++)
{
   $filename = $Files[$index];
   $chromosome = WhichChromAndSeq($filename);
   $PREFIXtoID = "$chromosome".'_';

   $n++;

   # Re-initialize for each chromosome
   $insertion_event_id = 0; # The number of 3' intact L1 insertions. 
   $AllInsertions = 0; # TOTAL number of L1 insertions.

   # In this next sub, put all relevant coordinates in an array
   # and merge Starts with Ends, if appropriate. Recall, due to 
   # the sometimes-insertion in the 5' UTR, the library used in
   # RepeatMasker has the starts and ends of the L1 split up. 
   # These are a special kind of merge since if there's a start
   # adjacent to an end and the distance is ~128 bp, allow them
   # to be merged.
   &MergeStEnd($filename);

   # From the file made in sub MergeStEnd, once again collect 
   # pertinent information and place in arrays.

   open(RMOUT2, "RMOut2") || die("Cannot open RMOut2 for reading: $!");

   $line_number = 0;
   while (defined($Line = <RMOUT2>))
   {
      $line_number++;
      ($gi[$line_number], $start[$line_number], $end[$line_number], $strand[$line_number], $Piece[$line_number], $L1start[$line_number], $L1end[$line_number], $L1left[$line_number], $PercId[$line_number], $ExtraSeq[$line_number]) = split (/\s+/, $Line);
   }
   close(RMOUT2) || die("Cannot close RMOut2: $!");

   # Calculate merged coordinates based on merge parameters, for each L1
   # Only report (after merging) L1s that end within 30 bp of end of 3' UTR
   # Create *.tab files for 5' and 3' flanking regions of defined lengths (see
   # defaults and user-defined parameters at the beginning)

   $inversionFLAG = '';
   
   # Loop through each segment in the StEnd Merge file
   # Check to see if they can be merged
   for ($count = 1; $count <= $line_number; $count++)
   {
      if ($inversionFLAG eq "INV")
      {
         $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
	 $inversionFLAG = '';
      }
      elsif ( defined($gi[$count + 1]) && ($gi[$count] eq $gi[$count + 1]) )
      {
         $deltagi = abs($start[$count + 1] - $end[$count]);
	 if ($deltagi < $DELTAGI)
         # L1 segments adjacent when mapped onto query  
	 {
#####################################################

	    if (($strand[$count] eq "+") && ($strand[$count + 1] eq "+"))
	    {
	    # Both segments go in the same direction (deletion, 
            # duplication, insertion, or enough point mutations 
            # to disrupt alignment with L1.3 at this location)
	       $deltaL1 = abs($L1start[$count + 1] - $L1end[$count]);
	       if ($deltaL1 < $DELTAL1)
	       {
	          #MERGE
                  # Get the lengths before changing them!
                  $Length1 = abs($start[$count] - $end[$count]) + 1; 
                  $Length2 = abs($start[$count + 1] - $end[$count + 1]) + 1;
                  $PercIden = ChangePercId($PercId[$count], $PercId[$count + 1], $Length1, $Length2);
                  $PercId[$count + 1] = $PercIden;
  		  $start[$count + 1] = $start[$count];
		  $L1start[$count + 1] = $L1start[$count];
                  $Piece[$count + 1] = "$Piece[$count]"." $Piece[$count + 1]";
                  $ExtraSeq[$count + 1] = "$ExtraSeq[$count]"." $ExtraSeq[$count + 1]";
	       }
	       else 
	       {
	          $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
	       }
	    }

#####################################################

	    elsif (($strand[$count] eq "C") && ($strand[$count + 1] eq "+"))
	    {
	    # The only class for once-5'-inverted(/deleted) elements
	       if ($L1end[$count] < $L1end[$count + 1]) 
               # 3'segment runs in the same direction as query
	       {
	          $deltaL1 = abs($L1start[$count + 1] - $L1end[$count]);
		  if ($deltaL1 < $DELTAL1)
		  {
		     print INVERSIONS "$gi[$count]\t$start[$count] .. $end[$count + 1]\tconsists of inverted L1.3 segments: $L1start[$count] .. $L1end[$count] and $L1start[$count + 1] .. $L1end[$count + 1]\n";
		     $inversionFLAG = 'INV';
		     #MERGE
                     # Get the lengths before changing them!
                     $Length1 = abs($start[$count] - $end[$count]) + 1; 
                     $Length2 = abs($start[$count + 1] - $end[$count + 1]) + 1;
                     $PercIden = ChangePercId($PercId[$count], $PercId[$count + 1], $Length1, $Length2);
                     $PercId[$count + 1] = $PercIden;
  		     $start[$count + 1] = $start[$count];
		     $L1start[$count + 1] = $L1start[$count];
                     $Piece[$count + 1] = "$Piece[$count]"." $Piece[$count + 1]";
                     $ExtraSeq[$count + 1] = "$ExtraSeq[$count]"." $ExtraSeq[$count + 1]";
		  }
		  else
		  {
		     $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
		  }
	       }
	       else  # 3' segment runs on C strand of GenBank record
	       {
	          $deltaL1 = abs($L1start[$count] - $L1end[$count + 1]);
		  if ($deltaL1 < $DELTAL1)
		  {
		     print INVERSIONS "$gi[$count]\t($start[$count] .. $end[$count + 1])\tconsists of inverted L1.3 segments: $L1start[$count + 1] .. $L1end[$count + 1] and $L1start[$count] .. $L1end[$count]\n";
		     $inversionFLAG = 'INV';
		     #MERGE
		     # Get the lengths before changing them!
                     $Length1 = abs($start[$count] - $end[$count]) + 1; 
                     $Length2 = abs($start[$count + 1] - $end[$count + 1]) + 1;
                     $PercIden = ChangePercId($PercId[$count], $PercId[$count + 1], $Length1, $Length2);
                     $PercId[$count + 1] = $PercIden;
  		     $start[$count + 1] = $start[$count];
                     $L1left[$count + 1] = $L1left[$count];
                     $Piece[$count + 1] = "$Piece[$count]"." $Piece[$count + 1]";
                     $ExtraSeq[$count + 1] = "$ExtraSeq[$count]"." $ExtraSeq[$count + 1]";
		     $L1end[$count + 1] = $L1end[$count];
		     $strand[$count + 1] = "C";	#for 3'ts direction
		  }
		  else
		  {
		     $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
		  }
	       }
	    }

#####################################################

	    elsif (($strand[$count] eq "C") && ($strand[$count + 1] eq "C"))
	    {
	       #both segments in the same orientation (see above)
	       $deltaL1 = abs($L1end[$count + 1] - $L1start[$count]);
	       if ($deltaL1 < $DELTAL1)
	       {
	          #MERGE
                  # Get the lengths before changing them!
                  $Length1 = abs($start[$count] - $end[$count]) + 1; 
                  $Length2 = abs($start[$count + 1] - $end[$count + 1]) + 1;
                  $PercIden = ChangePercId($PercId[$count], $PercId[$count + 1], $Length1, $Length2);
                  $PercId[$count + 1] = $PercIden;
  		  $start[$count + 1] = $start[$count];
		  $L1end[$count + 1] = $L1end[$count];
		  $L1left[$count + 1] = $L1left[$count];
                  $Piece[$count + 1] = "$Piece[$count]"." $Piece[$count + 1]";
                  $ExtraSeq[$count + 1] = "$ExtraSeq[$count]"." $ExtraSeq[$count + 1]";
	       }
	       else
	       {
	          $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
	       }
	    }  # End processing an inversion merge event.

#####################################################

	    else  # Two adjacent entries cannot be merged.
	    {
	       $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
	    }

         } # End merging 2 adjacent output lines representing a single
           # insertion event with indel <30
#####################################################
         else
         {
            $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);  # The indel is >30
         }

      }  # End working with entries from the same gi number

#####################################################

      else # Two adjacent entries in RM output that are neither INV nor of same gi. 
      {
         $KeepIt = ReportInsertionEvent($L1left[$count], $strand[$count], $start[$count], $end[$count], $gi[$count], $L1start[$count], $L1end[$count], $inversionFLAG, $PercId[$count], $Piece[$count], $ExtraSeq[$count]);
      }
      if ($KeepIt)
      {
         $Gi[$insertion_event_id] = $gi[$count];	 
         # These four are used to calculate coordinates for pre-insertion locus
         $Strand[$insertion_event_id] = $strand[$count];
         $Start_in_Gi[$insertion_event_id] = $start[$count];
         $End_in_Gi[$insertion_event_id] = $end[$count];
         $L1left[$count] =~ s/\(|\)//g;
         $L1left_after_merge[$insertion_event_id] = $L1left[$count];
         # For 3'ts evaluation
         $inversionFLAG[$insertion_event_id] = $inversionFLAG;
         $L1_length[$insertion_event_id] = $end[$count] - $start[$count] + 1;
      }
   }  # End going through all lines in a given *.out file

close MERGED;
close INVERSIONS;

$new_insertion_events = $insertion_event_id;
print "Total of $new_insertion_events 3' intact new insertion events\n";
print "Total of $AllInsertions insertion events\n";

# PART 2 of the program:  consider flanking sequence of each insertion
# and find TSDs.

$no_TSDs_counter = 0;
$standard_insertion_counter = 0;
$transduction_events_counter = 0;
$TSDs_wo_ATail = 0;
$PreInsFile = "PRE_INSERTION_LOCUSChr"."$chromosome";
open (OUTPUT1, ">$PreInsFile") || die("Cannot create $PreInsFile: $!");
$StatFile = "STATSChr"."$chromosome";
open (OUTPUT2, ">$StatFile") || die("Cannot create $StatFile: $!");

for ($ins_event = 1; $ins_event <= $insertion_event_id; $ins_event++)
{
   # Retrieve fasta sequences of L1 flanks to search for TSDs
   $a_file = $ins_event.'_5pr.tab';
   $Out_file = $ins_event.'_5pr.fa';
   &RangeToFasta($a_file, $Out_file);
   system("rm $a_file\n");

   $b_file = $ins_event.'_3pr.tab';
   $Out_file = $ins_event.'_3pr.fa';
   &RangeToFasta($b_file, $Out_file);
   system("rm $b_file\n");

   # Run blast2sequences
   $i_file = $ins_event.'_5pr.fa';
   $j_file = $ins_event.'_3pr.fa';
   $o_file = $ins_event.'.bl2seq';

   # Increase the word size to 9 to avoid error that too many
   # HSPs are being found and bl2seq is limited to 200 HSPs.
   system ("bl2seq -i $i_file -j $j_file -p blastn -g F -o $o_file -W 9 -F F -S 1 -d 3000 -e 1000.0") && die ("Problem running bl2seq: $!");
   system ("rm $i_file");  # No need to keep 5' flank of L1
	   
   # Evaluate bl2seq results
   $hspFLAG = 0;  # hspFLAG indicates in which phase program is in for collecting
   $HSPcounter = 0;     #information about the TSD
   $TSDscore = 0;
   $GotAsForEvent = 0; # Boolean to indicate whether putative PolyAs have been collected
   $PolyAScore = 0;
   $CurrentBestTSD = "NoTail"; # Always favor a TSD with a tail, but if there's not
                               # one available, take the best one lacking a tail (for
                               # standard insertions only).

   open (BL2SEQ, "$o_file") || die ("Cannot open $o_file: $!");
   while ($line = <BL2SEQ>)  # This loop searches for the best pair of HSPs, i.e. TSDs
   {
      if ($hspFLAG == 0)	
      {
         if ($line =~ /^\sIdentities = (\d+)\/(\d+)\s/)
	 {
	    $HSPcounter ++;
	    $HSPidentities = $1;
	    $HSPlength = $2;
	    $hspFLAG = 1 if ($HSPlength > 8);
	 }
	 elsif ($line =~ /^Lambda\s+/)  # end of data in bl2seq
	 {
	    print "$HSPcounter HSPs processed for insertion event number $ins_event\n";
	    if ($TSDscore == 0)	
	    {  # No HSP found with a polyA tail preceding the 3' match, OR no HSPs found.
	       $no_TSDs_counter ++;
	       print "No TSDs found for insertion event number $ins_event\n";
	    }
	    else  
            {  # At least 1 pair of HSPs scored YES for polyA tail -- Ver11: polyA not required
               if ($CurrentBestTSD eq "NoTail")
	       {
                  $TSDs_wo_ATail++ if ($CurrentBestTSD eq "NoTail");
	       }
               elsif ($transduction_state ne "3prTS")
	       {
	          $standard_insertion_counter++;
	       }
               elsif ($transduction_state eq "3prTS")
	       {
                  $transduction_events_counter++;
	       }

	       # Produce composite TSD sequence (case of mismatches between 5' and 3' HSPs)
	       if ($TSD_5pr_sequence eq $TSD_3pr_sequence)
	       {
	          $TSD = $TSD_5pr_sequence;
	       }
	       else
	       {
	          print "different TSD sequences: $TSD_5pr_sequence and  $TSD_3pr_sequence\t";
		  @Tsd_5pr = split (//, $TSD_5pr_sequence);
		  @Tsd_3pr = split (//, $TSD_3pr_sequence);
		  $tsd_length = @Tsd_5pr;
		  die("TSD length incompatibility in insertion event number $ins_event.  Die!\n") if ($tsd_length != @Tsd_3pr);
		  for ($bp = 0; $bp < $tsd_length; $bp++)
		  {
		  # print "$bp$Tsd_5pr[$bp]$Tsd_3pr[$bp].";
		     if ($Tsd_5pr[$bp] ne $Tsd_3pr[$bp])
		     {
		        print "$bp-th bp is different: $Tsd_5pr[$bp] vs. $Tsd_3pr[$bp]\t";
		        $Tsd_5pr[$bp] = '['.$Tsd_5pr[$bp].'/'.$Tsd_3pr[$bp].']';
		        print "new $bp-th character is $Tsd_5pr[$bp]\tComposite ";
		     }
		  }
	          $TSD = join ("", @Tsd_5pr);	
	       }
	       print "TSD sequence :  $TSD\n";
					
               # get gi and start/end coordinates for 
               # pre-insertion locus -- saved in sub 
               # ReportInsertionEvent retrieve flanking
               # sequences for pre-insertion locus output

	       open(preinsFIVEprFLANK, ">PREINS_5pr.tab") || die ("Cannot create PREINS_5pr.tab: $!");
	       open(preinsTHREEprFLANK, ">PREINS_3pr.tab") || die ("Cannot create PREINS_3pr.tab: $!");

	       if ($Strand[$ins_event] eq "+")
	       {
                  # 1st bp of 5'TSD, in the genbank record
	          $TSD_5pr_start_inGi = $Start_in_Gi[$ins_event] - $FIVE_PR_FLANK + $TSD_5pr_start - 1; 	 	
		  $preinsertion_locus_start_inGi = $TSD_5pr_start_inGi - $FLANKS_PREINSERTION;
                  # last bp of 3'TSD, in the genbank record
		  $TSD_3pr_end_inGi = $End_in_Gi[$ins_event] - $L1_TIPS + $TSD_3pr_end;    
		  $preinsertion_locus_end_inGi = $TSD_3pr_end_inGi + $FLANKS_PREINSERTION;
		  $pos_before_5prTSD = $TSD_5pr_start_inGi - 1;
		  $pos_after_3prTSD = $TSD_3pr_end_inGi + 1;
                  # These next two statements are for L1s that fall near the end
                  # of a sequence file.  They're not quite correct in the sequence
                  # they yield for the pre-insertion locus.
                  if ($preinsertion_locus_start_inGi < 0) 
		  {   $preinsertion_locus_start_inGi = 1;  }
                  if ($pos_before_5prTSD < 0)
		  {   $pos_before_5prTSD = 1;  }
		  print preinsFIVEprFLANK "$Gi[$ins_event]\t$preinsertion_locus_start_inGi .. $pos_before_5prTSD";
		  print preinsTHREEprFLANK "$Gi[$ins_event]\t$pos_after_3prTSD .. $preinsertion_locus_end_inGi";
		  close preinsFIVEprFLANK;
		  close preinsTHREEprFLANK;
                  &RangeToFasta("PREINS_5pr.tab", "PREINS_5pr.fa");
                  &RangeToFasta("PREINS_3pr.tab", "PREINS_3pr.fa");

		  open (FPR, "PREINS_5pr.fa") || die ("Cannot open PREINS_5pr.fa for insertion event $ins_event: $!");
		  @FIVE_PR = <FPR>;
		  close FPR;
		  open (TPR, "PREINS_3pr.fa") || die ("Cannot open PREINS_3pr.fa for insertion event $ins_event: $!");
		  @THREE_PR = <TPR>;
		  close TPR;
		  $preins_5pr_sequence = $FIVE_PR[1];
		  chomp $preins_5pr_sequence;
		  $preins_3pr_sequence = $THREE_PR[1];
		  chomp $preins_3pr_sequence;
                  system("rm PREINS_3pr.fa");
                  system("rm PREINS_5pr.fa");
                  # For output files, convert the id representing
                  # an intact 3' end to an id representing its 
                  # position wrt ALL L1 insertions.
                  $ConvertId = $AllVsIntact{$ins_event};
		  print OUTPUT1 ("$PREFIXtoID$ConvertId\t$Gi[$ins_event]\t$preinsertion_locus_start_inGi .. $preinsertion_locus_end_inGi\t$preins_5pr_sequence\t$TSD\t$preins_3pr_sequence\t$TSDscore\t$TSDPolyATail\t$TSDPolyATailMotif\t$pos_before_5prTSD\t$pos_after_3prTSD\t$transduction_state\n"); 
		  print OUTPUT2 ("$PREFIXtoID$ConvertId\t$preins_5pr_sequence\t$TSD\t$preins_3pr_sequence\t$TSDscore\t$transduction_state\t$L1_length[$ins_event]\t$inversionFLAG[$ins_event]\n");

 		  if ($transduction_state eq "3prTS")
		  {
		     $transduction_filename = $PREFIXtoID.$ConvertId.'_transd.tab';
		     open (TRANSDUCTION, ">$transduction_filename") || die ("Cannot create file $transduction_filename: $!");
		     $transd_segment_start = $End_in_Gi[$ins_event] + 1;
                     # Include 3' TSD in transduced segment sequence, as a
                     # control during further 3'ts sequence analysis
		     $transd_segment_end = $TSD_3pr_end_inGi;
                     # For ease in subsequent analysis, place both the 
                     # transduced sequence AND L1 coords in this file.	
		     print TRANSDUCTION ("$Gi[$ins_event]\t$transd_segment_start .. $transd_segment_end\n");
                     print TRANSDUCTION ("$Gi[$ins_event]\t$Start_in_Gi[$ins_event] .. $End_in_Gi[$ins_event]");
                     # There's a *_transduction.tab file for each transduced segment
		     close TRANSDUCTION; 
		  }
	       }

	       elsif ($Strand[$ins_event] eq "C")
	       {
	          $TSD_5pr_start_inGi = $End_in_Gi[$ins_event] + $FIVE_PR_FLANK - $TSD_5pr_start + 1;          # 1st bp of 5'TSD, in the genbank record
		  $TSD_3pr_end_inGi = $Start_in_Gi[$ins_event] + $L1_TIPS - $TSD_3pr_end;
                  # last bp of 3'TSD, in the genbank record
		  $preinsertion_locus_start_inGi = $TSD_3pr_end_inGi - $FLANKS_PREINSERTION;
		  $preinsertion_locus_end_inGi =  $TSD_5pr_start_inGi +  $FLANKS_PREINSERTION;
		  $pos_after_5prTSD = $TSD_5pr_start_inGi + 1; 
		  $pos_before_3prTSD = $TSD_3pr_end_inGi - 1; 
		  print preinsFIVEprFLANK "$Gi[$ins_event]\t($pos_after_5prTSD .. $preinsertion_locus_end_inGi)";
                  # These next two statements are for L1s that fall near the end
                  # of a sequence file.  They're not quite correct in the sequence
                  # they yield for the pre-insertion locus.
                  if ($preinsertion_locus_start_inGi < 0) 
		  {   $preinsertion_locus_start_inGi = 1;  }
                  if ($pos_before_3prTSD < 0)
		  {   $pos_before_3prTSD = 1;  }
		  print preinsTHREEprFLANK "$Gi[$ins_event]\t($preinsertion_locus_start_inGi .. $pos_before_3prTSD)";
		  close preinsFIVEprFLANK;    #.tab file
		  close preinsTHREEprFLANK;   #.tab file
                  &RangeToFasta("PREINS_5pr.tab", "PREINS_5pr.fa");
                  &RangeToFasta("PREINS_3pr.tab", "PREINS_3pr.fa");
		  open (FPR, "PREINS_5pr.fa") || die ("Cannot open PREINS_5pr.fa for insertion event $ins_event: $!");
		  @FIVE_PR = <FPR>;
		  close FPR;
		  open (TPR, "PREINS_3pr.fa") || die ("Cannot open PREINS_3pr.fa for insertion event $ins_event: $!");
		  @THREE_PR = <TPR>;
		  close TPR;
		  $preins_5pr_sequence = $FIVE_PR[1];
		  chomp $preins_5pr_sequence;
		  $preins_3pr_sequence = $THREE_PR[1];
		  chomp $preins_3pr_sequence;
                  system("rm PREINS_3pr.fa");
                  system("rm PREINS_5pr.fa");
                  # For output files, convert the id representing
                  # an intact 3' end to an id representing its 
                  # position wrt ALL L1 insertions.
                  $ConvertId = $AllVsIntact{$ins_event};
		  print OUTPUT1 "$PREFIXtoID$ConvertId\t$Gi[$ins_event]\t($preinsertion_locus_start_inGi .. $preinsertion_locus_end_inGi)\t$preins_5pr_sequence\t$TSD\t$preins_3pr_sequence\t$TSDscore\t$TSDPolyATail\t$TSDPolyATailMotif\t$pos_before_3prTSD\t$pos_after_5prTSD\t$transduction_state\n";
		  print OUTPUT2 "$PREFIXtoID$ConvertId\t$preins_5pr_sequence\t$TSD\t$preins_3pr_sequence\t$TSDscore\t$transduction_state\t$L1_length[$ins_event]\t$inversionFLAG[$ins_event]\n";

 		  if ($transduction_state eq "3prTS")
		  {
		     $transduction_filename = $PREFIXtoID.$ConvertId.'_transd.tab';
		     open (TRANSDUCTION, ">$transduction_filename") || die ("Cannot create file $transduction_filename: $!");
		     $transd_segment_start = $Start_in_Gi[$ins_event] - 1;
                     # Include 3' TSD in transduced segment sequence
                     # as a control during further 3'ts sequence analysis
		     $transd_segment_end = $TSD_3pr_end_inGi;
                     # For ease in subsequent analysis, place both the 
                     # transduced sequence AND L1 coords in this file.
		     print TRANSDUCTION ("$Gi[$ins_event]\t($transd_segment_end .. $transd_segment_start)\n");
                     print TRANSDUCTION ("$Gi[$ins_event]\t($Start_in_Gi[$ins_event] .. $End_in_Gi[$ins_event])");
		     close TRANSDUCTION;
		  }
	       }
					
	    }
	    last;
	 }
      }
      elsif ($hspFLAG == 1)  # Now need query (5'flank) positions for already found HSP
      {
         if ($line =~ /^Query: (\d+)\s+([a|t|g|c|n|-]+)\s+(\d+)/)
	 {
	    $HSP_5pr_start = $1;
	    $HSP_5pr_sequence = $2;
	    $HSP_5pr_end = $3;
	    $hspFLAG = 2;
	 }
      }
      else  # $hspFLAG = 2 -- Now need subject (3'flank) positions for the same HSP
      {
         if ($line =~ /^Sbjct: (\d+)\s+([a|t|g|c|n|-]+)\s+(\d+)/)
         {
	    $HSP_3pr_start = $1;
	    $HSP_3pr_sequence = $2;
	    $HSP_3pr_end = $3;

            if (!($GotAsForEvent))
	    {
               # Recall, $j_file is the 3' flanking sequence
               &FindPolyATails($j_file, $PREFIXtoID, $ins_event, $Gi[$ins_event]);
               $GotAsForEvent = 1; 
	    }

	    print "  $HSPcounter\t$HSPidentities of $HSPlength match\n";
	    print "   \t$HSP_5pr_start to $HSP_5pr_end in 5' flank\n";
	    print "   \t$HSP_3pr_start to $HSP_3pr_end in 3' flank\n";

            if ($Strand[$ins_event] eq "+")
	    { 
               $ThrPrTSDstart = $End_in_Gi[$ins_event] - $L1_TIPS + $HSP_3pr_start;
               $ThrPrTSDend = $End_in_Gi[$ins_event] - $L1_TIPS + $HSP_3pr_end; 
            }
            elsif ($Strand[$ins_event] eq "C")
            {  # These starts and ends are oriented wrt the + strand of the gi record. 
               # e.g. start is more 5' than the end      
               $ThrPrTSDstart = $Start_in_Gi[$ins_event] + $L1_TIPS - $HSP_3pr_end;
               $ThrPrTSDend = $Start_in_Gi[$ins_event] + $L1_TIPS - $HSP_3pr_start;
            }

	    @ReturnVals = AdjPolyATail($PREFIXtoID, $ins_event, $Gi[$ins_event], $ThrPrTSDstart, $ThrPrTSDend);
            $Result = $ReturnVals[0];
            $Score = $ReturnVals[1];
            $Motif = $ReturnVals[2];

            $TSDTail = $ReturnVals[3]; # The polyA tail affiliated with the TSD
            if ($Result eq "YES") # As in "Yes, the TSD is preceded by a polyA tail."
	    {
	       if ( ( (&TsdScore($HSP_5pr_end, $HSP_3pr_start, $HSPlength, $HSPidentities) >= $TSDscore) && (($Score >= $PolyAScore) || ($Motif)) ) || ($CurrentBestTSD eq "NoTail") )
	       {
                  if ( ($Motif) && ($Score < 10) )
		  {
                     # Mediocre effort to avoid low A bias in regular motifs. 
                     $Score = 10; 
		  }
                  $CurrentBestTSD = "YesTail";
                  $PolyAScore = $Score;
	          $TSDscore = &TsdScore($HSP_5pr_end, $HSP_3pr_start, $HSPlength, $HSPidentities);
	          $TSD_5pr_sequence = $HSP_5pr_sequence;  
	          $TSD_3pr_sequence = $HSP_3pr_sequence;  
	          $TSD_5pr_start = $HSP_5pr_start;	
	          $TSD_3pr_end = $HSP_3pr_end;
                  $TSDPolyATail = $TSDTail;
                  $TSDPolyATailMotif = $Motif;
                  $LengthOfTail = length($TSDPolyATail);
	          # Check for 3' transduction
	          if ($HSP_3pr_start > ($L1_TIPS + $L1left_after_merge[$ins_event] + 20 + $LengthOfTail) ) # at least 20 bp of unique DNA if this L1 broke off L1.3 within the last 30 bp AND had a comparable polyA tail after that
	          {
	             $transduction_state = '3prTS';
	          }
	          else
	          {
	             $transduction_state = '';
		  }
	       }
	    } # End 'if' loop for having found a valid polyA tail
            else  # No polyA tail was found 
	    {
               # Only allow standard insertions and do not replace a tail-containing TSD.
               # As a standard insertion, only allow 10 bp of wiggle room for the 
               # placement of the 3' TSD that lacks a polyA.
               if ( ($CurrentBestTSD ne "YesTail") && ($HSP_3pr_start < ($L1_TIPS + $L1left_after_merge[$ins_event] + 10)) && (&TsdScore($HSP_5pr_end, $HSP_3pr_start, $HSPlength, $HSPidentities) >= $TSDscore) )
	       {
                  $CurrentBestTSD = "NoTail";
                  $TSDscore = &TsdScore($HSP_5pr_end, $HSP_3pr_start, $HSPlength, $HSPidentities);
	          $TSD_5pr_sequence = $HSP_5pr_sequence;  
	          $TSD_3pr_sequence = $HSP_3pr_sequence;  
	          $TSD_5pr_start = $HSP_5pr_start;	
	          $TSD_3pr_end = $HSP_3pr_end;
                  $TSDPolyATail = "X";
                  $TSDPolyATailMotif = 0;
                  $LengthOfTail = length($TSDPolyATail);
                  $transduction_state = '';
               }
            }
	    $hspFLAG = 0;
	 }
      }
   } # End 'while' loop for looping through each blast2seq result
   system ("rm $j_file");  # cleanup
   system ("rm $o_file");  # cleanup
} # End 'for' loop for looping through each insertion event found by RM
 
# Remove the polyA files of previous chromosome so that they don't 
# grow and interfere with subsequent analyses.
system("rm PolyACoords PolyAMerge PolyAMergeFinal");
$MkDir = "TransdChr"."$chromosome";
system("mkdir $MkDir");
system("mv *_transd.tab $MkDir");

$MkMainDir = "Chr"."$chromosome";
system("mkdir $MkMainDir");
system("mv $MkDir $MkMainDir");
system("mv INVERSION* MERGE* STATS* PRE_INSERTION* $MkMainDir");

close OUTPUT1;
close OUTPUT2;
close TRANSDUCTION;

open (SUM, ">>SUMMARY1") || die ("Cannot create SUMMARY: $!");
print SUM "\nCHROMOSOME $chromosome\n";
print SUM "$AllInsertions L1 insertions were recorded.\n";
print SUM " $new_insertion_events are intact at 3' end\n";
print SUM "   $no_TSDs_counter of them have no TSDs\n";
print SUM "   $TSDs_wo_ATail of them have TSDs with poor polyA tail\n";
print SUM "   $standard_insertion_counter are standard insertion with A tail\n";
print SUM "   $transduction_events_counter are 3' transduction events\n";
close SUM;

&RemoveSeqFiles;

} # End going through all *.out files	

print "Processed $n RepeatMasker .out files \n";

} # End of main block of code
#############################################################

# SUBROUTINE DEFINITIONS FOLLOW

############################################################
sub FindPolyATails # At the sequence 3' of the L1.
{                  
   # A sequence will be searched through for putative polyA tails.   
   # Sequence in file is passed to subroutine. Same sequence 
   # file that was passed to Blast2Sequences as 3' flank.

   my $SeqToAnalyze = $_[0];              
   my $Prefix = $_[1];  
   my $InsertionNumber = $_[2];
   my $Gi = $_[3];

   # Continuous length of non-A nucs allowed in tail - "Contaminating island"
   my $CONTAMTHRESHOLD = 2; # Used to be 3
   # Lower bounds on percent As required to call a polyA tail short, long, etc.
   my $VShrtTHRESHOLD = 0.60;
   my $ShrtTHRESHOLD = 0.65;
   my $LngTHRESHOLD = 0.75; 
   my $VLngTHRESHOLD = 0.80;

   my $Sequence = ''; # Input sequence to be analyzed.

   # Obtain information about the incoming sequence
   my($Coord_line, $GiPiece, $FlankStart, $FlankEnd, $Strand);
   my($DNASeqLine, $Length, @Sequence);
   # Initialize some values
   my($StartedATail, $PolyASeq, $ContaminatingSeq, $MaxContam, $NumAs, $i);
   # PolyA tail characteristics
   my($Start, $PercAs, $RegularMotif, $Okay);
   # Store strings whose polyA status is still questionable
   my($MaybePolyA, $NewLength, $MaybeNumAs, $ChkAEffect, $AppendIt);
   # Special polyA motifs information
   my($MotifThreshold, $ALength, $ContamLength);

   open(DNASEQ, "$SeqToAnalyze") || die("Cannot open $SeqToAnalyze: $!");

   $Coord_line = <DNASEQ>; # First line should contain coordinates of sequence

   if ($Coord_line =~ /(\d+)\s+(\d+) .. (\d+)/)
   {
      $GiPiece = $1;
      $FlankStart = $2;
      $FlankEnd = $3;
      $Strand = "+";
   } 
   elsif ($Coord_line =~ /(\d+)\s+\((\d+) .. (\d+)\)/)
   {
      $GiPiece = $1;
      $FlankStart = $2;
      $FlankEnd = $3;
      $Strand = "C";
   }

   print "Coord_line='$Coord_line'\n";
   print "GiPiece=$GiPiece, FlankStart=$FlankStart, FlakEnd=$FlankEnd, Strand=$Strand\n";

   while ( !(eof(DNASEQ)) )
   {
      $DNASeqLine = <DNASEQ>;
      chomp($DNASeqLine);
      if ( !($DNASeqLine =~ /^>/) )  # Not the def line
      {
         $Sequence = $Sequence.$DNASeqLine;
      }   
   }

   close(DNASEQ) || die("Cannot close DNASEQ filehandle: $!"); 

   $Length = length($Sequence);
   @Sequence = split(//, $Sequence);  # Need input sequence as an array.

   $StartedATail = 0; # Boolean
   $PolyASeq = '';
   $ContaminatingSeq = '';
   $MaxContam = '';
   $NumAs = 0;
   $i = $L1_TIPS + 1; # Start at 11 b/c the first 10 bp in this 
                     # sequence are part of the 
                     # L1 3' UTR by virtue of the use of $L1TIPS = 10.

   open (POLYAOUT, ">PolyACoords") || die("Cannot open report file PolyACoords: $!");
   print POLYAOUT (">$Prefix$InsertionNumber\t$Gi\n");

   # As the sequence is traversed, can either start a polyA sequence, extend
   # a polyA sequence, collect contaminating nucs w/in a polyA sequence, stop
   # collecting nucs to add to a polyA sequence, or none of the above.

   while ($i < $Length)
   {
      if ( ($Sequence[$i] =~ /A/i) && (!($StartedATail)) )
      {
         $StartedATail = 1;
         $Start = $i + 1;  # B/c array starts at index 0
         $PolyASeq = $PolyASeq.$Sequence[$i];  
         $NumAs++;
         $PercAs = $NumAs/(length($PolyASeq) + length($ContaminatingSeq));
         $RegularMotif = 0; # Boolean defining e.g. AATAATAATAATAAT pattern
         $Okay = 1;
      }
      elsif ( ($Sequence[$i] =~ /A/i) && ($StartedATail) )
      {
         $MaybePolyA = $PolyASeq.$ContaminatingSeq.$Sequence[$i];
         $NewLength = length($MaybePolyA);
         $MaybeNumAs = $NumAs + 1;
         $ChkAEffect = $MaybeNumAs/length($MaybePolyA); # Effect on %As if added.
         $AppendIt = 0;
         if (length($ContaminatingSeq) == 0)# No contaminating sequence, add the A.
	 {
            $AppendIt = 1;
	 }
         elsif (($RegularMotif) && ($ChkAEffect >= $MotifThreshold))
	 {
            $AppendIt = 1 if (length($ContaminatingSeq) < $CONTAMTHRESHOLD);
	 }
         elsif (($NewLength < 5) && ($ChkAEffect > $VShrtTHRESHOLD))
         {
            ($AppendIt = 1) if (length($ContaminatingSeq) < $CONTAMTHRESHOLD);
	 }
         elsif (($NewLength < 10) && ($ChkAEffect > $ShrtTHRESHOLD))
         {
            ($AppendIt = 1) if (length($ContaminatingSeq) < $CONTAMTHRESHOLD);
	 }
         elsif (($NewLength < 30) && ($ChkAEffect >= $LngTHRESHOLD))
         {
            ($AppendIt = 1) if (length($ContaminatingSeq) < $CONTAMTHRESHOLD);
	 } 
         elsif (($NewLength >= 30) && ($ChkAEffect >= $VLngTHRESHOLD))
         {
            ($AppendIt = 1) if (length($ContaminatingSeq) < $CONTAMTHRESHOLD);
	 }
 
         if ($AppendIt)
         {
            $PolyASeq = $PolyASeq.$ContaminatingSeq.$Sequence[$i];
            $ContaminatingSeq = '';
            $NumAs++;
            $PercAs = $NumAs/(length($PolyASeq));
         }
         else
	 {
            $Okay = 0; # The current contam. sequence cannot be appended
            $i--;      # without crossing PercA threshold.  Abort current tail.
         }            
      }
      elsif ( (!($Sequence[$i] =~ /A/i)) && ($StartedATail) )   
      {
         $ContaminatingSeq = $ContaminatingSeq.$Sequence[$i];
         if ( ($PolyASeq =~ /(A+)([CTG])\1\2\1\2\1/) ) # Special motif pattern
	 {
            # Make sure it's not aborted due to low PercAs so tailor
            # the threshold with $MotifThreshold
            $RegularMotif = 1;  #Boolean
            $ALength = length($1);
            $ContamLength =  length($2);
            $MotifThreshold = $ALength/($ContamLength + $ALength);
         }
      }

      if ( ($StartedATail) && (!($Okay)) )  # Store the current found tail and re-set 
      { 
         if (length($PolyASeq) > 1) # Record tails of length 2 bp or longer
         {
            &Record($PolyASeq, $Strand, $NumAs, $FlankStart, $FlankEnd, $Start);
         }
         $StartedATail = 0;  # Re-set the variables to find another tail.
         $PolyASeq = '';
         $ContaminatingSeq = '';
         $MaxContam = '';
         $NumAs = 0;
      }
      $i++;
   }

   if ( ($StartedATail) && $Okay )  # EOF before storing a final tail
   { 
      if (length($PolyASeq) > 1) # Only record tails of length 3 bp or longer
      {
         &Record($PolyASeq, $Strand, $NumAs, $FlankStart, $FlankEnd, $Start);
      }
   }
   print POLYAOUT "//\n";
   close (POLYAOUT) || die("Cannot close PolyACoords at end of sub: $!");
   &MergeTails($Prefix, $InsertionNumber);
}
############################################################
sub Record
{
   my($PolyASeq, $Strand, $As, $FlankStart, $FlankEnd, $Start) = @_;
   print "Record('$PolyASeq', '$Strand', '$As', '$FlankStart', '$FlankEnd', '$Start')\n";

   my($Aend, $Astart, $TotalLength);

   if ($Strand eq "C")
   {
      $Aend = $FlankEnd - $Start + 1;
      $Astart = $Aend - length($PolyASeq) + 1;
      $TotalLength = length($PolyASeq);
      print POLYAOUT ("($Astart .. $Aend)  $PolyASeq  $As out of $TotalLength\n");
   }
   elsif ($Strand eq "+")
   {
      $Astart = $Start + $FlankStart - 1; # Convert to gi coordinates
      $Aend = $Astart + length($PolyASeq) - 1;
      $TotalLength = length($PolyASeq);
      print POLYAOUT ("$Astart .. $Aend  $PolyASeq  $As out of $TotalLength\n");
   }
}
###########################################################
sub MergeTails
{
   # This is a clean-up after having found putative polyA tails.  In this
   # sub, adjacent tails are tested to see if they can be merged into a 
   # long tail, providing the polyA content is at least 73%.
   # Also, tails that satisfy the Regular Motif ATATATAT (n.b. Just one 
   # A in between each contaminant) are pieced together if they fell
   # through the cracks above. 

   my $Prefix = $_[0];
   my $Insertion = $_[1];

   # Constants used
   my $MERGEPERCTHRESHOLD = 0.73; # Minimum percent As required for 2 tails to be merged.
   my $MINLENGTH = 10; # Minimum length of polyA tail to be reported.
   my $CONTAMTHRESHOLD = 2; # Used to be 3
   my $CONSIDEREDLONG = 20; # Merge tails if their combined length is at least this
                         # length and they're only 1 or 2 bases apart.

   # For reading in the PolyA tails found in sub FindPolyATails
   my($InLine, @WholeLine, @Starts, @Ends, @Tail, @NumAs, @Length, $Complement);
   my($ArrayLength, $i);

   # For deciding if two adjacent tails can be merged
   my($j, $Separation, $PutativeLength, $TotalAs, $Percent, $IncludingAMotif);
   my($ProbablyOkay, $NewRegMotif);

   open (POLYAFILE, "PolyACoords") || die("Cannot open PolyACoords: $!");
   open (NEWPOLYAFILE, ">PolyAMerge") || die("Cannot open PolyAMerge: $!");

   while (defined($InLine = <POLYAFILE>)) 
   {
      if ($InLine =~ /^>$Prefix$Insertion/)
      {
         print NEWPOLYAFILE ("$InLine");
         $i = 0;
         while (defined($InLine = <POLYAFILE>) && (!($InLine =~ /\/\//)))
         {
            if ($InLine =~ /(\d+) .. (\d+)  (\w+)  (\d+) out of (\d+)/) # + strand
	    {
               chomp($InLine);
               $WholeLine[$i] = $InLine;
               $Starts[$i] = $1;
               $Ends[$i] = $2;
               $Tail[$i] = $3;
               $NumAs[$i] = $4;
               $Length[$i] = $5;
               $Complement = 0;
            }
            elsif ($InLine =~ /\((\d+) .. (\d+)\)  (\w+)  (\d+) out of (\d+)/) # C strand
	    {
               chomp($InLine);
               $WholeLine[$i] = $InLine;
               $Starts[$i] = $1;
               $Ends[$i] = $2;
               $Tail[$i] = $3;
               $NumAs[$i] = $4;
               $Length[$i] = $5;
               $Complement = 1;
            }
            $ArrayLength = $i;
            $i++;  
         }
      }
   } # Finished filling arrays
   close (POLYAFILE) || die("Cannot close PolyACoords: $!");
   if ($Complement)
   {
      @WholeLine = reverse(@WholeLine);
      @Starts = reverse(@Starts);
      @Ends = reverse(@Ends);
      @Tail = reverse(@Tail);
      @NumAs = reverse(@NumAs);
      @Length = reverse(@Length);
   }

   for ($j = 0; $j < $ArrayLength; $j++)         
   {          
      $Separation = $Starts[$j + 1] - $Ends[$j] - 1;
      if ($Separation <= $CONTAMTHRESHOLD)  # Consider merging tails of up to 3 bp away.
      {
         $PutativeLength = $Length[$j] + $Length[$j + 1] + $Separation;
         $TotalAs = $NumAs[$j] + $NumAs[$j + 1];
         $Percent = $TotalAs/$PutativeLength;
         if ($Separation == 1)
	 {
            if ($IncludingAMotif) # Some motif tails are so dirty that they don't quite
	    {                     # fit the pattern below.  Thus, make some concessions
               $ProbablyOkay = 1; # for appending to them if the separation is a mere 
               $IncludingAMotif = 0; # 1 nucleotide.
            }
            if ( ($Tail[$j + 1] =~ /$Tail[$j]/) || ($Tail[$j] =~ /$Tail[$j + 1]/) )
	    {
               if ($NumAs[$j + 1] == ($Length[$j + 1] - 1))  # Only one non-A nucleotide.
	       {
                  $NewRegMotif = 1;
                  $IncludingAMotif = 1;
               }
            }
         }
         if ( ($Percent >= $MERGEPERCTHRESHOLD) || ($NewRegMotif) || ($ProbablyOkay) )
         { 
            $Starts[$j + 1] = $Starts[$j];
            $NumAs[$j + 1] = $TotalAs;
            $Length[$j + 1] = $PutativeLength;
            if ($Complement)
	    {
               $Tail[$j + 1] = $Tail[$j + 1]. "[$Separation]". $Tail[$j];
               $WholeLine[$j + 1] = "($Starts[$j] .. $Ends[$j + 1])  $Tail[$j + 1]  $TotalAs out of $PutativeLength";
	    }
            else
	    {
               $Tail[$j + 1] = $Tail[$j]. "[$Separation]". $Tail[$j + 1];
               $WholeLine[$j + 1] = "$Starts[$j] .. $Ends[$j + 1]  $Tail[$j + 1]  $TotalAs out of $PutativeLength";
            }
            $NewRegMotif = 0;
            $ProbablyOkay = 0;
	 }
         else
         {
            print NEWPOLYAFILE ("$WholeLine[$j]\n") if ($Length[$j] >= $MINLENGTH);
            $IncludingAMotif = 0;
         }
      }
      else
      {
         print NEWPOLYAFILE ("$WholeLine[$j]\n") if ($Length[$j] >= $MINLENGTH);
         $IncludingAMotif = 0;
      }
   }
   # Due to design of loop, the last entry is not included.
   print NEWPOLYAFILE ("$WholeLine[$ArrayLength]\n") if ($Length[$ArrayLength] >= $MINLENGTH);
   $IncludingAMotif = 0;

   close (NEWPOLYAFILE) || die("Cannot close PolyAMerge: $!");
   # Now, as a final merge, if two "long" tails are separated by just 1 or 2 
   # nucleotides, then merge them together. 

   # Re-use these array names, thus need to empty them.
   @WholeLine = ();
   @Starts = ();
   @Ends = ();
   @Tail = ();
   @NumAs = ();
   @Length = ();
   $ArrayLength = 0;

   open (NEWPOLYAFILE, "PolyAMerge") || die("Cannot open PolyAMerge: $!");
   open (FINALPOLYAFILE, ">PolyAMergeFinal") || die("Cannot open PolyAMergeFinal: $!");

   while (defined($InLine = <NEWPOLYAFILE>)) 
   {
      if ($InLine =~ /^>$Prefix$Insertion/)  # Beginning of an insert.
      {
         print FINALPOLYAFILE ("$InLine"); # Def'n line
         $i = 0;
         while (defined($InLine = <NEWPOLYAFILE>) && (!($InLine =~ /\/\//)))
         {
            chomp($InLine);
            $WholeLine[$i] = $InLine;
            if ($InLine =~ /\(\d+ .. \d+\)/) # Parentheses because complement
	    {
               $InLine =~ s/\(//;
               $InLine =~ s/\)//;  # Get rid of parentheses
               $Complement = 1;
            }
            else 
	    {
               $Complement = 0;
            }
            if ($InLine =~ /(\d+) .. (\d+)  (\S+)  (\d+) out of (\d+)/) 
	    {
               $Starts[$i] = $1;
               $Ends[$i] = $2;
               $Tail[$i] = $3;
               $NumAs[$i] = $4;
               $Length[$i] = $5;
            }
            $ArrayLength = $i;
            $i++;  
         }
      }
   } # Finished filling arrays
   close (NEWPOLYAFILE) || die("Cannot close PolyAMerge: $!");

   for ($j = 0; $j < $ArrayLength; $j++)         
   {          
      if (($j + 1) < $ArrayLength)
      {
         $Separation = $Starts[$j + 1] - $Ends[$j] - 1;
         $PutativeLength = $Length[$j] + $Length[$j + 1] + $Separation;
      }
      else
      {
         $PutativeLength = 0; # Prevent going thru the next merging loop
      }

      if ( ($Separation < $CONTAMTHRESHOLD) && ($PutativeLength >= $CONSIDEREDLONG) ) 
      {
         $TotalAs = $NumAs[$j] + $NumAs[$j + 1];
         $Starts[$j + 1] = $Starts[$j];
         $NumAs[$j + 1] = $TotalAs;
         $Length[$j + 1] = $PutativeLength;
         if ($Complement)
	 {
            $Tail[$j + 1] = $Tail[$j + 1]. "[$Separation]". $Tail[$j];
            $WholeLine[$j + 1] = "($Starts[$j] .. $Ends[$j + 1])  $Tail[$j + 1]  $TotalAs out of $PutativeLength";
         }
         else
	 {
            $Tail[$j + 1] = $Tail[$j]. "[$Separation]". $Tail[$j + 1];
            $WholeLine[$j + 1] = "$Starts[$j] .. $Ends[$j + 1]  $Tail[$j + 1]  $TotalAs out of $PutativeLength";
         }
      }
      else
      {
         print FINALPOLYAFILE ("$WholeLine[$j]\n");
      }

   }
   ## Due to design of loop, the last entry is not included.
   print "ArrayLength=$ArrayLength, WholeLine=$WholeLine[$ArrayLength]\n";
   print FINALPOLYAFILE ("$WholeLine[$ArrayLength]\n");
   close (FINALPOLYAFILE) || die("Cannot close FINALPOLYAFILE: $!");

}
           
############################################################
sub AdjPolyATail  # Check to see if the end of a putatative TSD is adjacent
{                 # to a putative polyA tail.
   my $Prefix = $_[0];
   my $Insertion = $_[1];
   my $Gi = $_[2];
   my $TSDstart = $_[3];
   my $TSDend = $_[4]; 

   my($MINLENGTH) = 10; # Minimum length of a polyA tail
   my($Input_line, $Strand, $Tail);
   my(@ScoreAndTail, $PolyAStart, $PolyAEnd); 
   my($Score, $TSDTail, $Motif);    
   my($NumAs, $Nuc, @Tail);

   $Result = "NO"; 

   open (POLYAMERGE, "PolyAMergeFinal") || die("Cannot open PolyAMergeFinal: $!");
   while ($Input_line = <POLYAMERGE>)
   {
      if ($Input_line =~ /$Prefix$Insertion\t$Gi/)
      {
         while (defined($Input_line = <POLYAMERGE>))
	 {
            if ($Input_line =~ /\((\d+) .. (\d+)\)  (\w+)/)
	    {
               $Strand = "C";
               $Input_line =~ s/\(//;
               $Input_line =~ s/\)//;
            }
            else
	    {
               $Strand = "+";
	    }
            if ($Input_line =~ /(\d+) .. (\d+)  (\w+)/)
	    {  
               $PolyAStart = $1;
               $PolyAEnd = $2; 
               $Tail = $3;
               # Next conditional statement requires that if TSD is embedded 
               # in a polyA tail, allow at least a length of 5 bp [what if 
               # it's 5 bp of junk?] for polyA tail alone OR if TSD is
               # positionned after the polyA, it can't be more than 2 bp away.

               if ($Strand eq "+")
	       {
                  if ( (($TSDstart - $PolyAStart >= $MINLENGTH) && ($TSDstart < $PolyAEnd)) || ((($TSDstart - $PolyAEnd) < 3) && ($TSDstart - $PolyAEnd) > 0) )
                  {
                     $Result = "YES";
                     # Now take the fraction of the tail that precedes 
                     # the TSD (if applicable) and report IT as the 
                     # PolyA tail and score it (score is number of As).
                     @ScoreAndTail = PolyAScore($PolyAStart, $PolyAEnd, $TSDstart, $TSDend, $Gi, $Strand); 
                     $Score = $ScoreAndTail[0];
                     $TSDTail = $ScoreAndTail[1];
                     @Tail = split(//, $TSDTail);
                     $NumAs = 0;
                     foreach $Nuc (@Tail)
                     {
                        if ($Nuc =~ /A/i)
                        {
                           $NumAs++;  # The total number of As will be the score.
                        }
                     }
                     $Result = "NO" if ($NumAs/@Tail < 0.75);
                     $Result = "NO" if (length($TSDTail) < 10);
                     $Motif = $ScoreAndTail[2];
                  } 
	       }
               elsif ($Strand eq "C")
               {
                  if ( (($PolyAEnd - $TSDend >= $MINLENGTH) && ($PolyAStart < $TSDend)) || ((($PolyAStart - $TSDend) < 3) && (($PolyAStart - $TSDend) > 0)) )
                  {
                     $Result = "YES";
                     # Now take the fraction of the tail that precedes 
                     # the TSD (if applicable) and report IT as the 
                     # PolyA tail and score it (score is number of As).
                     @ScoreAndTail = PolyAScore($PolyAStart, $PolyAEnd, $TSDstart, $TSDend, $Gi, $Strand); 
                     $Score = $ScoreAndTail[0];
                     $TSDTail = $ScoreAndTail[1];
                     @Tail = split(//, $TSDTail);
                     $NumAs = 0;
                     foreach $Nuc (@Tail)
                     {
                        if ($Nuc =~ /A/i)
                        {
                           $NumAs++;  # The total number of As will be the score.
                        }
                     }
                     $Result = "NO" if ($NumAs/@Tail < 0.73);
                     $Result = "NO" if (length($TSDTail) < 10);
                     $Motif = $ScoreAndTail[2];
                  } 
               }
            }
	 } # End going through all the candidate polyA tails
      }
      else
      {
         die("The PolyA tail file is off sync with the current insertion event.\n");
      }
   }
   close (POLYAMERGE);
   return ($Result, $Score, $Motif, $TSDTail);
}
############################################################
sub ReportInsertionEvent
{
   my($L1left, $strand, $start, $end, $gi, $L1start, $L1end, $inversionFLAG, $PercId, $Piece, $ExtraSeq) = @_;
   my($fpf_name, $tpf_name, $five_start, $five_end, $three_start);
   my($three_end, $ThrPrIntact); 
   my $KeepIt = 0; # Boolean indicating whether or not this 
                    # insertion event should be kept.
   $L1left =~ s/\(|\)//g;

   # Here, record ALL L1 insertions in the MERGE file -- not
   # just th 3' intact ones.
   if ($strand eq "+")
   {
      $ThrPrIntact = " ";
      if ( ($L1left <= $THR_PR_END) && ($Piece ne "FullStart") )
      {
         $KeepIt = 1;
         $ThrPrIntact = "ThrPrIntact";
         $insertion_event_id ++;
	 print "insertion_event_id = $insertion_event_id\n";
	 $fpf_name = $insertion_event_id.'_5pr.tab';
	 open (FIVEprFLANK, ">$fpf_name") || die("Cannot create $fpf_name: $!");
	 print "Created file $fpf_name\n";
	 $tpf_name = $insertion_event_id.'_3pr.tab';
	 open (THREEprFLANK, ">$tpf_name") || die("Cannot create $tpf_name: $!");
	 print "Created file $tpf_name\n";

	 $five_start = $start - $FIVE_PR_FLANK;
	 if ($five_start < 0) 
            { $five_start = 1; } # Avoid problems with negative nos. in RangeToFasta
	 $five_end = $start + $L1_TIPS - 1;
	 print FIVEprFLANK "$gi\t$five_start .. $five_end";
	 $three_start = $end - $L1_TIPS + 1;
	 $three_end = $end + $THREE_PR_FLANK;
	 print THREEprFLANK "$gi\t$three_start .. $three_end";
      }
      $AllInsertions++;
      print MERGED "$PREFIXtoID$AllInsertions\t$gi\t$start .. $end\t$L1start\t$L1end\t$L1left\t$PercId\t$ThrPrIntact\t$Piece\t$ExtraSeq\t$inversionFLAG\n";
      if ($ThrPrIntact eq "ThrPrIntact")
      {
         $AllVsIntact{$insertion_event_id} = $AllInsertions;
      }
   }

   elsif ($strand eq "C")
   {
      $ThrPrIntact = " ";
      if ( ($L1left <= $THR_PR_END) && ($Piece ne "FullStart") )
      {
         $KeepIt = 1;
         $ThrPrIntact = "ThrPrIntact";
         $insertion_event_id ++;
	 print "insertion_event_id = $insertion_event_id\n";
	 $fpf_name = $insertion_event_id.'_5pr.tab';
	 open (FIVEprFLANK, ">$fpf_name") || die ("Cannot create $fpf_name: $!");
	 print "Created file $fpf_name\n";
	 $tpf_name = $insertion_event_id.'_3pr.tab';
	 open (THREEprFLANK, ">$tpf_name") || die ("Cannot create $tpf_name: $!");
	 print "Created file $tpf_name\n";

	 $five_start = $end - $L1_TIPS + 1;
	 $five_end = $end + $FIVE_PR_FLANK;
	 print FIVEprFLANK "$gi\t($five_start .. $five_end)";
	 $three_start = $start - $THREE_PR_FLANK;
	 if ($three_start < 0) 
	    { $three_start = 1; } # Avoid problems w/ negative nos. in RangeToFasta
	 $three_end = $start + $L1_TIPS - 1;
	 print THREEprFLANK "$gi\t($three_start .. $three_end)";
      }      
      $AllInsertions++;
      print MERGED "$PREFIXtoID$AllInsertions\t$gi\t($start .. $end)\t$L1start\t$L1end\t$L1left\t$PercId\t$ThrPrIntact\t$Piece\t$ExtraSeq\t$inversionFLAG\n";
      if ($ThrPrIntact eq "ThrPrIntact")
      {
         $AllVsIntact{$insertion_event_id} = $AllInsertions;
      }
   }
   close FIVEprFLANK;
   close THREEprFLANK;

   return($KeepIt);
}

###########################################################################

sub TsdScore
{
   my $HSP_5pr_end = $_[0];
   my $HSP_3pr_start = $_[1];
   my $HSPlength = $_[2];
   my $HSPidentities = $_[3];

   my $return_score = 0;
   my $fpf_position_score = 0;
   my $tpf_position_score = 0;
   my $hsp_quality_score = 0;
   my $acceptFLAG = '';
   my($max_mismatches, $mismatches);
	
   if ($HSP_5pr_end > ($FIVE_PR_FLANK - $L1_TIPS))  # assign 5' HSP position score
   {
      $fpf_position_score = 100;
   }
   else  # linear score function with max = 80
   {
      $fpf_position_score = $HSP_5pr_end * 80 / ($FIVE_PR_FLANK - $L1_TIPS) ;
   }
	
   # Admittedly, this + 40 is somewhat random b/c don't know the length
   # of the polyA tail.  
   if ($HSP_3pr_start < ($L1_TIPS + 40) )   # assign 3' HSP position score 
   {
      $tpf_position_score = 20;
   }

   if ($HSPlength < 11)      # 8, 9, 10
   {       
      if ($HSPlength == $HSPidentities)
      {   
         $hsp_quality_score = 10;
      }
      else  # mismatch in short HSP
      {
         $acceptFLAG = 'NO';
      }
   }
   elsif (($HSPlength >= 11) && ($HSPlength <= 18))
   {
      if ($HSPlength == $HSPidentities)
      {
         $hsp_quality_score = 20;
      }
      elsif ( ($HSPlength - $HSPidentities) == 1)
      {
         $hsp_quality_score = 10;
      }
      else  # more than one mismatch in 11-18 bp HSPs
      {
         $acceptFLAG = 'NO';
      }
   }
   else # HSP length > 18 bp
   {
      $max_mismatches = int ($HSPlength / 50) + 1;
      $mismatches = $HSPlength - $HSPidentities;
      if ($mismatches <= $max_mismatches)
      {
         $hsp_quality_score = 10;
      }
      else  # more than allowed number of mismatches (1/50 bp)
      {
         $acceptFLAG = 'NO';
      }
   }

   $return_score = $fpf_position_score + $tpf_position_score + $hsp_quality_score  unless ($acceptFLAG eq 'NO');

   $return_score;
}

##################################################################

# This sub splits the contig sequence file into sub-files:  one file
# for each gi.  This should be faster than searching this master 
# sequence file for each gi each time.  Call these files the gi number
# with a .seq suffix.  Need to include the header line for subsequent
# BLAST searches.  Need to chomp off the newline characters to avoid
# counting them when collecting a subsequence.

sub SplitSeqFiles
{
   my($ContigFile, $DNASeqLine, @LineItems, $FoundGiNumber);

   while ($ContigFile = <*.fa>)  # This suffix represents the contig files
   {                               # filled with gi records in fasta format.
      open(DNASEQ, "$ContigFile") || die("Cannot open $ContigFile: $!");
      $DNASeqLine = <DNASEQ>;
      while ( !(eof(DNASEQ)) )
      {
         if ($DNASeqLine =~ /^>gi/)
         {
            @LineItems = split(/\|/, $DNASeqLine);
            $FoundGiNumber = $LineItems[1];
            open(NEWGIFILE, ">$FoundGiNumber".".seq");
            print NEWGIFILE ("$DNASeqLine\n");
            chomp($DNASeqLine = <DNASEQ>);
            while ((!($DNASeqLine =~ /^>gi/)) && (!(eof(DNASEQ))) )
            {
               print NEWGIFILE ("$DNASeqLine");
               chomp($DNASeqLine = <DNASEQ>);
            }
         }
         close(NEWGIFILE);
      }
      system("rm $ContigFile"); # Save space
      close(DNASEQ);    
   }
}
##################################################################

# Clean up.  Remove all the gi sequence files.

sub RemoveSeqFiles

{
   system("rm *.seq");
}

##################################################################

# Input is a gi number and a start nucleotide and an end nucleotide.
# If this range of nucleotides is in parentheses, this represents a
# sequence on the complement strand.  This sub returns the DNA sequence
# from 5' to 3' end that this range represents.

sub RangeToFasta
{
   my $FileName = $_[0];
   my $OutFile = $_[1];

   my $Complement = 0;
   my $Header = '';
   my($SeqRangeInfo, $GiOfInterest, $StartNucl, $EndNucl);
   my($RangeOfNucs, $Skip, $SeqFile, $Input_line);
   my($DefLine, @SeqArray, $SeqOfInterest, $LinesOfSeq, $i);


   open(RANGEINFO, "$FileName") || die("Cannot open $FileName: $!");
   if (defined(chomp($SeqRangeInfo = <RANGEINFO>)))
   {
      if ($SeqRangeInfo =~ /(\d+)\s+(\d+)\s\.\.\s(\d+)/) 
      {   
         $GiOfInterest = $1;
         $StartNucl = $2;
         $EndNucl = $3;
      }
      elsif ($SeqRangeInfo =~ /(\d+)\s+\((\d+)\s\.\.\s(\d+)\)/) 
      {   
         $GiOfInterest = $1;
         $StartNucl = $2;
         $EndNucl = $3;
         $Complement = 1;
      }
   }
   close(RANGEINFO) || die("Cannot close $FileName: $!");

   open(OUTPUT, ">$OutFile") || die("Cannot open $OutFile: $!");

   $RangeOfNucs = $EndNucl - $StartNucl + 1;
   $Skip = $StartNucl - 1; # Number of nucs to skip over before starting
                           # the collection.
   $SeqFile = $GiOfInterest.".seq";
   open(DNASEQ, "$SeqFile") || die("Cannot open $SeqFile: $!");
   
   if (($Input_line = <DNASEQ>) =~ /^>gi/)
   {
      chomp($Input_line);
      $Header = $Input_line;  # The definition line.
   }
   else
   {
      print("There's something wrong with the header of file $SeqFile.\n");
   }

   seek(DNASEQ, $Skip, 1);  # Place pointer at the next line which is DNA sequence.
   read(DNASEQ, $SeqOfInterest, $RangeOfNucs);  # Collects the sub-sequence of interest.

   if ($Complement)
   {
      $SeqOfInterest = MakeComplement($SeqOfInterest);
   }

   $DefLine = ">".$SeqRangeInfo;
   $DefLine =~ s/\t/ /g;
   print OUTPUT ("$DefLine\n");
   @SeqArray = $SeqOfInterest =~ /\w{1,60}/g; # Break up sequence into lines 
   $LinesOfSeq = @SeqArray;                   # of 1-60 nucleotides.
   for ($i = 0; $i < $LinesOfSeq; $i++)
   { print OUTPUT ("$SeqArray[$i]\n");}

   close(OUTPUT) || die("Cannot close $OutFile: $!");
   close(DNASEQ) || die("Cannot close $SeqFile: $!");

} # End sub

#################################################################
sub MakeComplement
{
   my $Sequence = $_[0];

   my(@SeqAsArray, $Length, $index, $NewNucl);
   my $NewSequence = "";

   @SeqAsArray = split(//, $Sequence);  # Make the sequence an array of chars.
   $Length = @SeqAsArray;
   for ($index = ($Length - 1); $index >= 0; $index--)
   {
      if ($SeqAsArray[$index] =~ /A/i)  # Ignore case.
      {  
         $NewNucl = "T";  
         $NewSequence = $NewSequence.$NewNucl;
      }
      elsif ($SeqAsArray[$index] =~ /T/i)
      {  
         $NewNucl = "A";  
         $NewSequence = $NewSequence.$NewNucl;
      }
      elsif ($SeqAsArray[$index] =~ /C/i)
      {
         $NewNucl = "G";  
         $NewSequence = $NewSequence.$NewNucl;
      }
      elsif ($SeqAsArray[$index] =~ /G/i)
      {  
         $NewNucl = "C";  
         $NewSequence = $NewSequence.$NewNucl;
      }
      elsif ($SeqAsArray[$index] =~ /N/i)
      {  
         $NewNucl = "N";  
         $NewSequence = $NewSequence.$NewNucl;
      }
   }
   return $NewSequence;
} # End sub    
############################################################
sub PolyAScore
{
   # This sub will fetch the polyA tail associated with a putative
   # TSD (n.b. Sometimes a fraction of the polyA tail immediately
   # upstream of a found HSP will be deemed part of the TSD.).  The
   # score is simply the total number of As in the tail.

   my($polyAstart, $polyAend, $tsdStart, $tsdEnd, $gi, $strand) = @_;

   my $TheTail = "";
   my($TheorPolyAStart, $TheorPolyAEnd, $TheorPolyA, $InLine, $Length);
   my(@Tail, $i, $MotifTag, $NumAs, $TempString, $CutOff);

   open(ASTRETCH, ">ARange") || die("Cannot open ASTRETCH filehandle: $!");
   if ($strand eq "C")
   {
      $TheorPolyAStart = $tsdEnd + 1;
      print ASTRETCH ("$gi\t($TheorPolyAStart .. $polyAend)");
   }
   elsif ($strand eq "+")
   {  
      $TheorPolyAEnd = $tsdStart - 1;
      print ASTRETCH ("$gi\t$polyAstart .. $TheorPolyAEnd");
   }   
   close(ASTRETCH) || die("Cannot close ARange: $!");
   &RangeToFasta("ARange", "ATail");
   open(ATAIL, "ATail") || die("Cannot open ATail in PolyAScore sub: $!");
   while (defined($InLine = <ATAIL>))
   {
      chomp($InLine);
      if (!($InLine =~ /^>/)) # Not the header line
      {
         $TheTail = $TheTail.$InLine;
      }
   }
   close(ATAIL) || die("Cannot close ATail in PolyAScore sub: $!");
   $Length = length($TheTail);
   @Tail = split(//, $TheTail);
   @Tail = reverse(@Tail);
   $NumAs = 0;
   $TempString = "";
   $CutOff = 0;
   for ($i = 0; $i < $Length; $i++)
   {
      # Go through tail backwards to counts the As and while doing
      # so, check to make sure that any clipping of the tail due to
      # TSD donation didn't leave it ending in non-As.  Since no more
      # than 2 contaminants in a row are allowed, only have to check 
      # the first 2 bases ($i == 0 and $i == 1).
      if ( ($Tail[$i] =~ /A/i) )
      {
         $NumAs++;  # The total number of As will be the score.
         $TempString .= $Tail[$i];
      }
      elsif ( ($Tail[$i] =~ /[CTG]/i) && ($i == 0) )
      {
         $CutOff = 1; 
      }
      elsif ( ($Tail[$i] =~ /[CTG]/i) && ($i == 1) && ($CutOff) )
      {
         # Do nothing
      }
      else
      {
         $TempString .= $Tail[$i];
      }
   }
   $TheTail = reverse($TempString);
   
   system("rm ARange ATail");
   if ($TheTail =~ /(A+)([CTG])\1\2\1\2\1/)
   {
      $MotifTag = 1;  
   }
   else
   {
      $MotifTag = 0;
   }
   return($NumAs, $TheTail, $MotifTag);
}

##############################################################
sub MergeStEnd
{   
   my($OutFile) = @_;

   my($line_number, $Line, $count, $i);
   my($skip0, @score, @PercSubs, @PercDels, @PercIns, @long_name, @start, @end);
   my($skip1, @strand, @Piece, $skip2, @coord1, @L1end, @coord3, @star, @PercIdentity);
   my(@getgi, @gi, @L1start, @L1left);
   my($SameGi, $SameStrand, $CorrectOrder); 
   my($Length1, $Length2, $NewPercId, $ExtraSeq);

   open(RMOUT, "$OutFile") || die("Cannot open $OutFile: $!");
   
   $line_number = 0;
   # Collect all pertinent info from RepeatMasker .out file into arrays.
   while (defined($Line = <RMOUT>))
   {
      $line_number++;
      ($skip0, $score[$line_number], $PercSubs[$line_number], $PercDels[$line_number], $PercIns[$line_number], $long_name[$line_number], $start[$line_number], $end[$line_number], $skip1, $strand[$line_number], $Piece[$line_number], $skip2, $coord1[$line_number], $L1end[$line_number], $coord3[$line_number], $star[$line_number]) = split (/\s+/, $Line);
      # Consider the PercentSubstitution field in RM output.  
      $PercIdentity[$line_number] = 100 - ($PercSubs[$line_number]);
   }
   close(RMOUT) || die("Cannot close $OutFile: $!");

   print "$line_number segments recognized in $OutFile\n";

   # The following arrays can now be used:
   # @score keeps scores from each line
   # @long_name keeps the name of sequence annotated (gi + extra symbols, here)
   # @start keeps start coordinate within genbank record (always lower than end)
   # @end keeps end coordinate within genbank record
   # @strand keeps C/+ strand info for L1 orientation in the genbank record
   # @coord1 keeps first L1 annotation coordinate
   # @L1end keeps second L1 annotation coordinate
   # @coord3 keeps third L1 annotation coordinate
   # @star keeps * information (whether a better L1 segment overlaps this one)

   # Extract fields needed for further scoring and sequence retrieval
   # @gi keeps gi for each line
   # @L1start keeps the start coordinate for each L1 (before 2nd Merge)
   # @L1left keeps the number (with parentheses) of L1 3' leftover nucleotides

   $count = 1;

   open(RMOUT2, ">RMOut2") || die("Cannot open RMOut2: $!");
   while ($count <= $line_number)
   {
      @getgi = split (/\|/, $long_name[$count]);
      $gi[$count] = $getgi[1];

      # Reverse orientation for ease in fetching 
      # a range of seqence from gi records
      if ($strand[$count] eq "C")
      {
       	 $L1start[$count] = $coord3[$count];
         $L1left[$count] = $coord1[$count];
      }
      elsif ($strand[$count] eq "+")
      {
      	 $L1start[$count] = $coord1[$count];
	 $L1left[$count] = $coord3[$count];
      }
      else {print "No meaningful strand annotation!!\n"};

      $count ++;
   }

   for ($i = 1; $i <= $line_number; $i++)
   {
      print RMOUT2 ("$gi[$i]   $start[$i]   $end[$i]   $strand[$i]   $Piece[$i]   $L1start[$i]   $L1end[$i]   ($L1left[$i])   $PercIdentity[$i]\n");
   } 
   close(RMOUT2) || die("Cannot close RMOut2: $!");
}
##################################################
sub ChangePercId
{
   # If merging two pieces, need to re-figure percent identity
   # to reference sequence since length changes.

   my($PercId1, $PercId2, $Length1, $Length2) = @_;
   my($FullLength);
 
   $FullLength = $Length1 + $Length2;

   $NewPerc = ((($PercId1/100) * $Length1) + (($PercId2/100) * $Length2)) / $FullLength;
   $NewPerc = $NewPerc * 100;
   return($NewPerc);
}

#################################################
sub WhichChromAndSeq
{
   my($filename) = @_;

   my($chromosome, $ChromDir, $MergeFile, $InversionFile, $SeqFile);

   if ($filename =~ /hs_chr(\S{1,2})\.fa\.out/)
   {
      $chromosome = $1;
   }
   else
   {
      print("Cannot recognize the chromosome from the filename.\n");
      print("Enter chromosome number:\n");
      $chromosome = <STDIN>;
      chomp $chromosome;
   }

   if ($chromosome =~ /\d{2}/)  # Two digit chromosome number
   {
      $ChromDir = "CHR_".$chromosome;
   }
   elsif ($chromosome =~ /\d{1}/)  # One digit chromosome
   {
      $ChromDir = "CHR_"."0".$chromosome;
   }
   else  # A non-digit - named chromosome
   {
      $ChromDir = "CHR_".$chromosome;
   } 

   print ("\nWorking on chromosome $chromosome\n");
   $MergeFile = "MERGEChr"."$chromosome";
   open (MERGED, ">$MergeFile") || die("Cannot create $MergeFile: $!");

   $InversionFile = "INVERSIONChr"."$chromosome";
   open (INVERSIONS, ">$InversionFile") || die("Cannot create $InversionFile: $!");

   $SeqFile = $filename; # The sequence filename is almost identical to the RM filename
   $SeqFile =~ s/\.out//;
   $SeqFile .= ".gz";

   # Go to the appropriate dir and fetch the chrom. build of interest

   system("cp $DirOfGenomeSeqs/$SeqFile $SeqFile");
   system("gunzip $SeqFile");
   &SplitSeqFiles;  # Take the 'master' contig sequence used for the Repeat
                 # Masker phase of analysis and split it up into files
                 # each containing a single DNA gi sequence.

   return($chromosome);
}

###################################################
sub IntroductoryRemarks
{

   my($Change_def, $number_of_arguments);

   print("The default values for this program currently are:\n");
   print("5' prime flank length for searching for TSDs = $FIVE_PR_FLANK bp\n");
   print("3' prime flank length for searching for TSDs = $THREE_PR_FLANK bp\n");
   print("Length of flanking sequence for preinsertion target reconstruction = $FLANKS_PREINSERTION bp\n");

   print("Do you want to alter these defaults? Y or N\n");
   $Change_def = <STDIN>;
   if ($Change_def =~ /Y/i)
   {
      print("To alter the default values, re-start this program.\n");
      print("On the command line, enter the revised numbers.  For example:\n");
      print("RepeatMasker   5' prime flank length for searching for TSDs   3' flank length for searching for TSDs  length of preinsertion flanking sequence \n");
      print("Enter all or none.\n");
      die;
   }

   $number_of_arguments = @ARGV;

   if ($number_of_arguments == 3)
   {
       print("$number_of_arguments arguments, reassigning the defaults..\n");
       $FIVE_PR_FLANK = $ARGV[0];
       $THREE_PR_FLANK = $ARGV[1];
       $FLANKS_PREINSERTION = $ARGV[2];
   }
   elsif (($number_of_arguments != 0) && ($number_of_arguments != 3))
   {
       print("Please enter the correct number of arguments! \n");
       print("Length of 5' flank (bp) for TSD finding, length of 3' flank for \n");
       print("TSD finding (bp), length of flanks to retrieve for pre-insertion locus\n");
       die;
   }

   print("Program will use the following values:\n");
   print("FIVE_PR_FLANK = $FIVE_PR_FLANK\n");
   print("THREE_PR_FLANK = $THREE_PR_FLANK\n");
   print("FLANKS_PREINSERTION = $FLANKS_PREINSERTION\n");
}	

