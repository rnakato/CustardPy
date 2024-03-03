#!/usr/bin/awk -f
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Takes a SAM file, looks for chimeric reads and normal reads.  Outputs
# only those reads mapped to chromosomes 1-24 with MAPQ > 0.
#
# Output to fname1 is of the form:
# strand1 chromosome1 position1 frag1 strand2 chromosome2 position2 frag2
#  mapq1 cigar1 seq1 mapq2 cigar2 seq2
#   where chr1 is always <= chr2
#
# Output to fname2 retains SAM format.
#
# Chimeric reads are treated as follows:
# "Normal" chimeras (with 3 reads, two with the ligation junction), where the
# third read (w/o ligation junction) is within 20Kbp of one of the ligation
# junction reads, are sent to fname1.  All other chimeras (where either there
# are more than 3 reads, or they don't contain the ligation junction, or the
# one non-ligation junction end is far from the others, are sent to fname2.
#
# awk -f chimeric_blacklist.awk -v fname1="norm_chimera" fname2="abnorm_chimera" fname3="unmapped"
# Juicer version 1.5

# returns absolute value
function abs(value)
{
    return (value<0?-value:value);
}
# returns minimum of two values
function min(value1,value2)
{
    return (value1<value2?value1:value2);
}
# examines read1 (s1,c1,p1) versus read2 and returns true if
# the first read comes before the second read
# this is so duplicates can be found after sorting even when the strand and
# chromosome are the same
function less_than(s1,c1,p1,s2,c2,p2)
{

    if (c1 < c2) return 1;
    if (c1 > c2) return 0;
    # c1 == c2
    if (s1 < s2) return 1;
    if (s1 > s2) return 0;
    # s1 == s2 && c1 == c2
    if (p1 < p2) return 1;
    if (p1 > p2) return 0;
    # all are equal, doesn't matter
    return 1;
}

function print_data(str, chr, pos, m, cigarstr, seq, name, fname, index1, index2) {
    if (less_than(str[index1], chr[index1], pos[index1], str[index2], chr[index2], pos[index2])) {
        print str[index1], chr[index1], pos[index1], str[index2], chr[index2], pos[index2], m[index1], cigarstr[index1], seq[index1], m[index2], cigarstr[index2], seq[index2], name[index1], name[index2] > fname;
    } else {
        print str[index2], chr[index2], pos[index2], str[index1], chr[index1], pos[index1], m[index2], cigarstr[index2], seq[index2], m[index1], cigarstr[index1], seq[index1], name[index2], name[index1] > fname;
    }
}

function calculate_adjusted_length(cigar) {
    seqlength = 0;
    currstr = cigar;
    # Count all matches, not just the first, to handle cases like 15M10S20M
    where = match(currstr, /[0-9]+[M|D|N|X|=]/);
    while (where > 0) {
        seqlength += substr(currstr, where, RLENGTH - 1) + 0;
        currstr = substr(currstr, where + RLENGTH);
        where = match(currstr, /[0-9]+[M|D|N|X|=]/);
    }
    # Add soft clipped bases to the position for proper adjustment
    if (cigar ~ /[0-9]+S$/) {
        where = match(cigar, /[0-9]+S$/);
        cigloc = substr(cigar, where, RLENGTH - 1) + 0;
        seqlength += cigloc;
    }
    # Optionally, add hard clipped bases as well
    else if (cigar ~ /[0-9]+H$/) {
        where = match(cigar, /[0-9]+H$/);
        cigloc = substr(cigar, where, RLENGTH - 1) + 0;
        seqlength += cigloc;
    }
    return seqlength;
}

# Function to adjust position based on strand and CIGAR string
function adjust_position(strand, cigar, pos) {
    split(cigar, cigar_parts, /[SH]/); # Split CIGAR string on S or H
    leading_bases = cigar_parts[1] + 0; # Convert leading bases count to number

    if (strand == 0 && cigar ~ /^[0-9]+S/) {
        pos -= leading_bases;
        if (pos <= 0) {
            pos = 1;
        }
    }
    else if (strand == 0 && cigar ~ /^[0-9]+H/) {
        pos -= leading_bases;
        if (pos <= 0) {
            pos = 1;
        }
    }
    else if (strand == 16) {
        # Call previously defined function to calculate adjusted length from CIGAR
        adjusted_length = calculate_adjusted_length(cigar);
        pos += adjusted_length - 1;
    }
    return pos; # Return the adjusted position
}

# Function to extract and assign data from a split array
function extract_and_assign_data(record, str, chr, pos, m, cigarstr, seq, name, j) {
    split(record, tmp); # Split the record into an array
    str[j] = and(tmp[2], 16); # Extract strand information using bitwise AND
    chr[j] = tmp[3]; # Chromosome
    pos[j] = tmp[4]; # Position
    m[j] = tmp[5]; # Mapping quality
    cigarstr[j] = tmp[6]; # CIGAR string
    seq[j] = tmp[10]; # Sequence
    # qual[j] could be extracted here if needed: qual[j] = tmp[11];
    name[j] = tmp[1]; # Read name
}

function calculate_distance(chr, pos, i, j) {
    return abs(chr[i] - chr[j]) * 10000000 + abs(pos[i] - pos[j]);
}

function check_mapped_read1and2(mapped, str, chr, pos, m, cigarstr, seq, name, read1, read2, header, c, count_norm, count_unmapped, fname1, fname3) {
	if (mapped[read1] && mapped[read2]) {
		count_norm++;
		print_data(str, chr, pos, m, cigarstr, seq, name, fname1, read1, read2);
	} else {
		if (count_unmapped == -1) {
			print header > fname3;
		}
		for (i in c) {
			print c[i] > fname3;
		}
		count_unmapped++;
	}
}

# Function to analyze chromatin interaction data and classify read pairs
function analyze_chromatin_interactions(chr, pos, mapped, c, header, fname1, fname2, fname3, count_unmapped, count_abnorm, count_norm) {
    # Calculate distances
    for (i=1; i<=4; i++) {
        for (j=i+1; j<=4; j++) {
            dist[i*10 + j] = calculate_distance(chr, pos, i, j)
        }
    }
    
    # Determine interaction pattern
    if ((dist[13] < 1000 && dist[24] < 1000) || (dist[14] < 1000 && dist[23] < 1000)) {
        read1 = 1; read2 = 2;
    } else if (dist[12] < 1000 && dist[34] < 1000) {
        read1 = 1; read2 = 3;
    } else {
        read1 = 0;
    }
    
    # Process based on determined pattern
    if (read1 != 0) {
		check_mapped_read1and2(mapped, str, chr, pos, m, cigarstr, seq, name, read1, read2, header, c, count_norm, count_unmapped, fname1, fname3) 
    } else {
        # Chimeric read with the 4 ends > 1KB apart
        if (count_abnorm == -1) {
            print header > fname2;
        }
        for (i in c) {
            print c[i] > fname2;
        }
        count_abnorm++;
    }
}

# Function to analyze distance metrics and classify read pairs
function classify_reads_based_on_distance(chr, pos, read, mapped, header, fname1, fname3, fname2, count_unmapped, count_abnorm, count_norm) {
    dist[12] = calculate_distance(chr, pos, 1, 2);
    dist[23] = calculate_distance(chr, pos, 2, 3);
    dist[13] = calculate_distance(chr, pos, 1, 3);
#    print dist[12], dist[23], dist[13]

    # Evaluate if any pair is within threshold
    if (min(dist[12], min(dist[23], dist[13])) < 1000) {
        # Determine read pairs based on specific conditions
        if (read[1] == read[2]) {
            read2 = 3;
            read1 = dist[13] > dist[23] ? 1 : 2;
        } else if (read[1] == read[3]) {
            read2 = 2;
            read1 = dist[12] > dist[23] ? 1 : 3;
        } else if (read[2] == read[3]) {
            read2 = 1;
            read1 = dist[12] > dist[13] ? 2 : 3;
        } else {
            printf("reads strange\n") > "/dev/stderr";
            exit 1;
        }
#        print read1, read2

        # Process based on mapping status
		check_mapped_read1and2(mapped, str, chr, pos, m, cigarstr, seq, name, read1, read2, header, c, count_norm, count_unmapped, fname1, fname3) 
    } else {
        # Process as chimeric read
        if (count_abnorm == -1) {
            print header > fname2;
        }
        count_abnorm++;
        for (i in c) {
            print c[i] > fname2;
        }
    }
}

# Function to extract read information and determine pair status
function extract_read_info(c_record, read, j) {
    split(c_record, tmp);
    split(tmp[1], readname, "/");
    # Check for backward compatibility and determine if read ends match
    if (length(readname) > 1) {
        read[j] = readname[2];
    } else {
        # Determine if it is the first in pair based on flag 64
        read[j] = (and(tmp[2], 64) > 0);
    }
}

BEGIN{
    OFS="\t";
    tottot = -1; # will count first non-group
}
$0 ~ /^@/{
    # print SAM header to SAM files
    header = header""$0"\n";
}
$0 !~ /^@/{
    if (tottot == -1) {
		header = header"@PG\tID:Juicer\tVN:1.6";
    }

    # input file is sorted by read name.  Look at read name to group
    # appropriately
    split($1,a,"/");

    print $0
#    print $1, a[1];
    if (a[1]==prev) {
		# move on to next record.  look below this block for actions that occur
		# regardless.
		count++;
#        print count
    } else {
		# deal with read pair group
		tottot++;
		if (count==3 || count==4) {
			# chimeric read
			for (j=1; j <= count; j++) {
				extract_read_info(c[j], read, j)
				extract_and_assign_data(c[j], str, chr, pos, m, cigarstr, seq, name, j);

               print "before", pos[j], chr[j], str[j]
				pos[j] = adjust_position(str[j], tmp[6], pos[j]);
                print "after", pos[j], chr[j], str[j]
				# blacklist - if 3rd bit set (=4) it means unmapped
				mapped[j] = and(tmp[2],4) == 0;
			}
			if (count == 4) {
				analyze_chromatin_interactions(chr, pos, mapped, c, header, fname1, fname2, fname3, count_unmapped, count_abnorm, count_norm)
			} else {
				classify_reads_based_on_distance(chr, pos, read, mapped, header, fname1, fname3, fname2, count_unmapped, count_abnorm, count_norm)
			}
		} else if (count > 3) {
			# chimeric read > 3, too many to deal with
			if (count_abnorm == -1) {
				print header > fname2;
			}
			for (i in c) {
				print c[i] > fname2;
			}
			count_abnorm++;
		} else if (count == 2) {
			# code here should be same as above, but it's a "normal" read
			j = 0;
			for (i in c) {
				extract_and_assign_data(c[i], str, chr, pos, m, cigarstr, seq, name, j);

				# blacklist - if 3rd bit set (=4) it means unmapped
				mapped[j] = and(tmp[2],4) == 0;

				pos[j] = adjust_position(str[j], tmp[6], pos[j]);
				j++;
			}
			if (mapped[0] && mapped[1]) {
				count_reg++;
				print_data(str, chr, pos, m, cigarstr, seq, name, fname1, 0, 1)
			}
			else {
				if (count_unmapped == -1) {
					print header > fname3;
				}
				for (i in c) {
					print c[i] > fname3;
				}
				count_unmapped++;
			}
		} else if (count == 1) {
			# this actually shouldn't happen, but it happens with alternate aligners on occasion
			if (count_abnorm == -1) {
				print header > fname2;
			}
			count_abnorm++;
			for (i in c) {
				print c[i] > fname2;
			}
		}
		# reset variables
		delete c;
		count=1;
		prev=a[1];
    }
    # these happen no matter what, after the above processing
    c[count] = $0;
}

END{
    # deal with read pair group
    tottot++;
    if (count==3 || count==4) {
		# chimeric read
		for (j=1; j <= count; j++) {
			extract_read_info(c[j], read, j)
			extract_and_assign_data(c[j], str, chr, pos, m, cigarstr, seq, name, j);

			pos[j] = adjust_position(str[j], tmp[6], pos[j]);
			# blacklist - if 3rd bit set (=4) it means unmapped
			mapped[j] = and(tmp[2],4) == 0;

		}
		if (count == 4) {
			analyze_chromatin_interactions(chr, pos, mapped, c, header, fname1, fname2, fname3, count_unmapped, count_abnorm, count_norm)
		} else {
			classify_reads_based_on_distance(chr, pos, read, mapped, header, fname1, fname3, fname2, count_unmapped, count_abnorm, count_norm)
		}
    } else if (count > 3) {
		# chimeric read > 3, too many to deal with
		if (count_abnorm == -1) {
			print header > fname2;
		}
		count_abnorm++;
		for (i in c) {
			print c[i] > fname2;
		}
    } else if (count == 2) {
		# code here should be same as above, but it's a "normal" read
		j = 0;
		for (i in c) {
			extract_and_assign_data(c[i], str, chr, pos, m, cigarstr, seq, name, j);

			# blacklist - if 3rd bit set (=4) it means unmapped
			mapped[j] = and(tmp[2],4) == 0;

			pos[j] = adjust_position(str[j], tmp[6], pos[j]);
			j++;
		}
		if (mapped[0] && mapped[1]) {
			count_reg++;
			print_data(str, chr, pos, m, cigarstr, seq, name, fname1, 0, 1)
		} else {
			if (count_unmapped == -1) {
			print header > fname3;
			}
			for (i in c) {
			print c[i] > fname3;
			}
			count_unmapped++;
		}
    }
    else if (count == 1) {
		# this actually shouldn't happen, but it happens with alternate aligners on occasion
		if (count_abnorm == -1) {
			print header > fname2;
		}
		count_abnorm++;
		for (i in c) {
			print c[i] > fname2;
		}
    }
    resfile=fname1".res.txt";
    printf("%d %d %d %d %d\n", tottot, count_unmapped, count_reg, count_norm, count_abnorm) >> resfile;
}
