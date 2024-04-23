#!/usr/bin/perl
use strict;
use warnings;

# 					To Do's						#
# 1. Allow commandline entry of the input file	#
# 2. Add simple logging							#
# 3. Add verbose logging						#
# 4. Add logging cleanup function				#
# 5. Format output (numbers mostly)				#
# ############################################# #
 
#                Initialize values                 #
# ################################################ #
my $inputfile = 'final_set_of_exons.fasta' ;
my $outputfile = 'final_set_of_exons_formatted.fasta.2' ;

my $DNAOutStringLength = 80;			# This is the max. number of characters each line of DNA will have.  HybPiper likes 80.
my $DNAString = "" ;					# The entire DNA string for each gene or gene segment
my $LastReadHeader = "" ;				# Saved the DNA header we just read from the input file, it will be the "next" gene we will process 
my $CurrentHeader = "" ;				# The DNA header of the gene we currently processing
my $DNAHeaderInCnt = 0 ;
my $DNAHeaderOutCnt = 0 ;
my $ReadCnt = 0 ;
my $DNAInCnt = 0 ;
my $DNAOutCnt = 1 ;
my $ValidDNA = "[ACGTURYSWKMBDHVN ]" ;
my $InValidDNACnt = 0 ;

# Open the input file in read mode #
open(my $in, '<', $inputfile) or die "Cannot open input file: $!";

# Open the output file in write mode #
open(my $out, '>', $outputfile) or die "Cannot open output file: $!";

# Open the Invalid Gene log file in write mode #
open(my $InvDNAlogout, '>', "InvalidDNA.log" ) or die "Cannot open output file: $!";


# Read from the input file and write formatted records to the output file #
while (my $line = <$in>) {
	$ReadCnt += 1 ;
	
#	print "0. ", $ReadCnt, "-Line read is: ", $line, "\n" ;
#	my $pause = <STDIN> ;
	# ######################################################### #
	# The first and last DNAHeader records need special processing	#
	
	if ($line =~ /^>/) {       # Is this a header record?
		$LastReadHeader = $line; 
		$DNAHeaderInCnt += 1;
				
		# 					For troubleshooting					#
		# ##################################################### #
#		print "1. LastReadHeader: ", $LastReadHeader, "\n";
#		print "2. CurrentHeader: ", $CurrentHeader, "\n" ;
#			$pause = <STDIN> ;
		
		# Is this the first header we are going to process?
		if ( $CurrentHeader eq "" ) {
			# We are reading our very first header, just need to save it and continue processing
			$CurrentHeader = $LastReadHeader ;
#			print "3. CurrentHdr = LastReadHdr: ", $CurrentHeader, " ", $LastReadHeader, "\n" ;
#			$pause = <STDIN> ;
		} elsif ($CurrentHeader ne $LastReadHeader) {
			
			# We have found a new header. We need to process the data we have already collected	#
			# #################################################################################	#
			
					# 					For troubleshooting					#
					# ##################################################### #
#				print "4. CurrentHdr ne LastReadHdr, Process the gene", "\n" ;
#				print "  CurrentHeader: ", $CurrentHeader, "\n", 
#					  "  LastReadHeader: ", $LastReadHeader, "\n",
#					  "  DNAOutStringLength: ", $DNAOutStringLength, "\n",
#					  "	 DNAString: ", $DNAString, "\n" ;
#				$pause = <STDIN> ;

				($CurrentHeader, $DNAString) = ScrubDNA( $CurrentHeader, $DNAString );
				$InValidDNACnt += IsDNAValid( $ValidDNA, $CurrentHeader, $DNAString ) ;
				WriteFormattedData( $CurrentHeader, $DNAString, $DNAOutStringLength ) ;
				
				# Set up for processing the gene	#
				# #################################	#
				$CurrentHeader = $LastReadHeader ;
				$DNAString = "" ;
		}
	}
	else {
			
		# This record is a DNA record not a head.				 				#
		# If there are multiple DNA records we will concatenate all of them		#
		# ##################################################################### #
		
				# 					For troubleshooting					#
				# ##################################################### #	
#		print "5. Saving DNA \n" ;

		$DNAString .= $line ;		# concatenate the lines of DNA for this gene
		$DNAInCnt += 1;
	}
}

# Processing of the last set of DNA		#
# #####################################	#

# 					For troubleshooting					#
# ##################################################### #
# Process the last set of DNA
#print "99. Processing last set of DNA", "\n" ;
#print "CurrentHdr ne LastReadHdr, Process the gene", "\n" ;
#print " CurrentHeader: ", $CurrentHeader, "\n", 
#	" LastReadHeader: ", $LastReadHeader, "\n",
#	" DNAOutStringLength: ", $DNAOutStringLength, "\n",
#	" DNAString: ", $DNAString, "\n" ;
# my $pause = <STDIN> ;

WriteFormattedData( $CurrentHeader, $DNAString, $DNAOutStringLength ) ;

# Close the filehandles
close($in);			# the FASTA file we processed
close($out);		# the new formatted FASTA file

print $InvDNAlogout "\nGenes with invalid DNA: ", $InValidDNACnt, "\n";
close($InvDNAlogout);	# the Invalid DNA log

print "#######################################################################\n" ;
print "  FASTA file to format: $inputfile\n";
print "  Formatted FASTA file: $outputfile\n";
print "    Total records read: ", $ReadCnt, "\n" ;
print "      DNA Headers read: ", $DNAHeaderInCnt, "\n" ;
print "              DNA read: ", $DNAInCnt, "\n" ;
print " Total records written: ", $DNAHeaderOutCnt + $DNAOutCnt, "\n" ;
print "   DNA Headers written: ", $DNAHeaderOutCnt, "\n" ;
print "           DNA written: ", $DNAOutCnt, "\n" ;
print "Genes with invalid DNA: ", $InValidDNACnt, " ** (see InvalidDNA.log)\n";
print "#######################################################################\n\n" ;
exit(0);


sub ScrubDNA {
		# Do some clean up of the DNA data we saved from the input file	#
		# #############################################################	#
	my( $myHeader, $myDNA ) = @_ ;

		my $SpaceCnt = 0 ;
			
		# Let's get ride of any line break characters #
		$myDNA =~ s/\r\n?|\n//g ;

		# Let's get ride of any spaces too. #
		if ($myDNA =~ s/\s+//g) {
			$SpaceCnt = 1 ;
		} ;
		print "Found a space in the DNA" if $SpaceCnt == 1;

		# Fix Header	#
	$myHeader = substr($myHeader, 0, rindex($myHeader, "-")) . "\n";  # strips the length value in the header
	$myHeader =~ s|(.+)-|$1.|;	# replaces the last occurance of "-" with a "."

return $myHeader, $myDNA  ;
}

sub IsDNAValid {
	
	my( $myValidDNA, $myHeader, $myDNA ) = @_ ;
	my $InValidDNAFound = 0 ;
	my $IfPrinted = 0 ;
	
	if( $myDNA !~ /^$myValidDNA+$/ ) {
    	if ($IfPrinted == 0 ) {
    		print $InvDNAlogout "Characters other than ", $myValidDNA, "found in DNA string.\n\n";
    		$IfPrinted = 1 ;
    }
    print $InvDNAlogout "Header is: ", $myHeader ;
    print $InvDNAlogout "   DNA is: ", $myDNA, "\n\n" ;
    	
    $InValidDNAFound = 1 ;
   }
   
return( $InValidDNAFound ) ;	
}

sub WriteFormattedData { 
		# The subroutine will take care of writing the header and all of the DNA data we have saved	#
		# from the input file to the formatted output file											#
		# #########################################################################################	#
	
	my( $Header, $DNA, $DNALineLength ) = @_ ;

		# 					For troubleshooting					#
		# ##################################################### #
#		print "W1. DNAHeader is: ", $Header, "\n" ;
#		print "W2. DNA is: ", $DNA, "\n" ;
#		print "W3. DNALineLength is: ", $DNALineLength, "\n" ;	
		my $FromHere = 0 ;


		my $DNALength = length($DNA) ;
#		print "DNALength: ", $DNALength, "\n" ;
#		my $ConsoleHeader = $Header;
#		$ConsoleHeader =~ s/\r\n?|\n//g ;
#		print "\rSaving: ", $ConsoleHeader ;		
		print $out $Header ;
		$DNAHeaderOutCnt += 1 ;
		
		my $DNAOutString = "" ;		# Reset it for the data we about to process
	
		# Loop throught the save DNA data writing it to the formatted output file	#
		# using a standard record length.											#
		# #########################################################################	#  										
		while ( $DNALength >= $FromHere ) {
			# 					For troubleshooting					#
			# ##################################################### #
#			print "FromHere: ", $FromHere, "\n" ;
#			my $pause2 = <STDIN> ;
			
			$DNAOutString = substr( $DNA, $FromHere, $DNALineLength) . "\n";

#			print "DNAOutString: ", $DNAOutString, "\n" ;
			
			print $out $DNAOutString ;
			$DNAOutCnt += 1 ;
			# 					For troubleshooting					#
			# ##################################################### #
#			print $out "\n" ;
#			print "W4. DNAOutString is: ", $DNAOutString, "\n" ;

			$FromHere += $DNAOutStringLength ;					# Update the starting position for getting the next set of DNA
#			print "W5. FromHere now is: ", $FromHere, "\n" ;
		}
		
		# 					For troubleshooting					#
		# ##################################################### #
#		my $pause = <STDIN> ;
		
	return ;
}