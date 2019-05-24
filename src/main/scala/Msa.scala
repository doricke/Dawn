
import scala.collection.mutable.{Map, ArrayBuffer, Set}

//  Author:     Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
//  Copyright:  Copyright (c) 2015 Massachusetts Institute of Technology
//  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

class Msa( val fileName: String, val name: String, wordSize: Int ) {

  val align = Map[String, Map[Int, Char]]()		// MSA

  val common_mers = Map[String, Set[String]]()		// [name, Set[k-mers]]

  val conservation = Map[Int, Float]()			// C(N) where +N = # of genes; N.M where M is number of classes; -N where N is number of differences

  val consensus = Map[String, Map[Int, Char]]()		// Consensus sequences

  val full_align = Map[String, Map[Int, Char]]()	// MSA - insertions added and gaps aligned.

  val gaps = Map[String, Map[Int, Gap]]()		// gap start-end

  val gene_names = Map[String, Int]()			// [gene names, count]

  val inserts = Map[String, Map[Int, Insert]]() 	// after:insertion

  val loops = Map[Int, Loop]()				// Group gaps and insert by loop region

  val mer_counts = Map[String, Int]()			// [k-mer, count]

  val observed = Map[String, Map[Int, Map[Char, Int]]]()	// [gene, [residue, count]]

  val proteins = Map[String, Protein]()			// [name, Protein]

  val rules = Array(Set('I', 'L', 'V'), Set('R', 'K'), Set('S', 'T'), Set('S', 'A'), Set('E', 'D'), Set('D', 'N'), Set('Q', 'E'), Set('Q', 'N'), Set('F', 'Y'))
  // val rules = Array(Set('I', 'L', 'V', '.'), Set('R', 'K', '.'), Set('S', 'T', '.'), Set('S', 'A', '.'), Set('E', 'D', '.'), Set('D', 'N', '.'), Set('Q', 'E', '.'), Set('Q', 'N', '.'), Set('F', 'Y', '.'))

  val segments = Map[String, Map[Int, Insert]]() 	// after:segment
 
  val seq_tools = new SeqTools()

  val tallies = Map[String, Map[Int, Char]]()		// Pre-alignment tallies to identify guide template

  val template = new StringBuilder()			// Consensus sequence identity template

  val top_mers = Set[String]()				// Set k-mers

  // ***************************************************************************
  def determineConservation() {
    common_mers foreach {case(seq_name, mers) =>
      if ( align contains seq_name ) {
        val gene = if ( proteins( seq_name ).fasta.annotation contains "gene" ) proteins( seq_name ).fasta.annotation( "gene" ).toUpperCase() else "other"
        val taxonomy = if ( proteins( seq_name ).fasta.annotation contains "taxonomy" ) proteins( seq_name ).fasta.annotation( "taxonomy" ).toUpperCase() else "Unknown"

        // if ( taxonomy.contains( "MAMMALIA" ) ) {
          if ( ( gene_names contains gene ) == false )
            gene_names += gene -> 1
          else
            gene_names( gene ) += 1
  
          // println( seq_name + " gene: " + gene )
          if ( ( observed contains gene ) == false )
            observed += gene -> Map[Int, Map[Char, Int]]()
   
          // println( "**** determine: " + seq_name + " size: " + align(seq_name).size + "\t" + align(seq_name).mkString )
   
          // for ( i <- 1 until (align( seq_name ).size+1) ) {
          for ( i <- 1 until (align( name ).size+1) ) {
            if ( ( observed( gene ) contains i ) == false )
              observed( gene ) += i -> Map[Char, Int]()
    
            // Add this residue to the count for this gene
            if ( align( seq_name ) contains i ) {
              val residue = align( seq_name )( i ).toUpper
              if ( residue != 'X' ) {
                if ( observed( gene )( i ) contains residue ) 
                  observed( gene )( i )( residue ) += 1
                else
                  observed( gene )( i ) += residue -> 1
              }
            }  // if
          }  // for
        // }  // if
        // else
          // println( "*** gene: " + gene + ", taxonomy: " + taxonomy )
      }  // if
    }  // foreach    

    computeConservation
  }  // determineConservation

  // ***************************************************************************
  def conservative( observed: Map[Char, Int] ): Char = {
    if ( observed.keys.size == 0 ) 
      '.'
    else {
      var cons = '.'
      rules foreach {case(rule) =>
        val diff_set = observed.keySet -- rule
        if ( diff_set.size == 0 ) {
          cons = observed.toList.sortBy{_._2}.last._1.toLower
        }  // if
      }  // foreach
      cons
    }  // if
  }  // conservative

  // ***************************************************************************
  def maxObserved( counts: Map[Char, Int] ): Tuple2[Char, Int] = {
    counts.toList.sortBy{_._2}.last
  }  // maxObserved

  // ***************************************************************************
  def computeConservation() {
    val con_all = Map[Int, Map[Char, Int]]()		// [position, [residue, count]]

    observed foreach {case(gene, tallies) =>
      consensus += gene -> Map[Int, Char]()

      // for ( i <- 1 until (tallies.keys.size+1) ) {
      for ( i <- 1 until (align(name).size+1) ) {
        if ( ( con_all contains i ) == false )  
          con_all += i -> Map[Char, Int]()

        if ( tallies contains i )
          if ( tallies( i ).keys.size == 1 )  {
            val residue = tallies( i ).keySet.head
            consensus( gene ) += i -> tallies( i ).keySet.head

            if ( residue != '.' )
              if ( ( con_all( i ) contains residue ) == false )
                con_all( i ) += residue -> 1
              else
                con_all( i )( residue ) += 1
          }
          else {
            val residue = conservative( tallies( i ) ).toUpper
            // consensus( gene ) += i -> conservative( tallies( i ) )
            consensus( gene ) += i -> residue

            if ( residue != '.' )
              if ( ( con_all( i ) contains residue ) == false )
                con_all( i ) += residue -> 1
              else
                con_all( i )( residue ) += 1
          }  // else
        else
          consensus( gene ) += i -> '.'
      }  // for
    }  // foreach

    println
    println( "All consensus" )
    for ( i <- 1 until (align(name).size+1) ) {
      val (residue, count) = if ( con_all( i ).keys.size > 0 ) con_all( i ).toList.sortBy{_._2}.last else new Tuple2('.', 0)
      // println( i + "  " + con_all( i ).mkString( "|" ) )
      val consensus_limit = if ( gene_names.size > 20 ) gene_names.size/4 else 5
      if ( count >= consensus_limit ) {
        template += residue
        // println( i + "\t" + residue + "\t" + count )
      }
      else {
        template += ' ' 
        // println( i + "\t\t\t" + residue + "\t" + count )
      }  // if
    }  // for
    println( "Consensus\t" + template + "\tConsensus" )
    val core_blocks = new CoreBlocks( template.toString )
    core_blocks.makeBlocks
    // optimizeCore( core_blocks.blocks )

    // Compute the conservation for the target query gene.
    gene_names foreach {case(gene_name, count) =>
      print( gene_name + "\t" )
      for ( i <- 1 until (consensus( gene_name ).keys.size+1) ) 
        print( consensus( gene_name )( i ) )
      println( "\t" + gene_name )
    }  // foreach

    // println
    // println( "Conservation level" )
    val gene = if ( proteins( name ).fasta.annotation contains "gene" ) proteins( name ).fasta.annotation( "gene" ).toUpperCase() else "other"
    for ( i <- 1 until (consensus( gene ).keys.size+1) ) {
      val res = consensus( gene )( i ).toUpper
      val res_set = Set( res )
      if ( res != '.' ) {
        conservation += i -> 1
        consensus foreach {case(gene_name, gene_cons) =>
          if ( ( gene != gene_name ) && ( consensus( gene_name ) contains i ) ) {
            val residue = consensus( gene_name )( i )
            if ( res == residue )
              conservation( i ) += 1
            else {
              val residue_set = Set( residue )
              rules foreach {case(rule) =>
                if ( ( (rule & res_set ).size > 0 ) && ( ( rule & residue_set ).size > 0 ) ) {
                  conservation( i ) += 1
                  // println( "*** match: " + res + ":" + residue )
                }
              }  // foreach
            }  // if
          }  // if
        }  // for
        // println( gene + "\t" + i + "\t" + res + "\t" + conservation( i ).toInt.toString )
      } else {
        // val variability = if ( observed(gene)(i) contains '-' ) 0 else -(observed( gene )( i ).size - 1)
        val variability = if ( observed(gene)(i) contains '-' ) -(observed(gene)(i).size) else -(observed( gene )( i ).size - 1)
        conservation += i -> variability
        // println( gene + "\t" + i + "\t" + res + "\t" + variability + "\t" + observed( gene )( i ).mkString( "|" ) )
      }  // if
    }  // for
    // println
    // println
  }  // computeConservation

  // ***************************************************************************
  def tallyMers() {
    common_mers foreach {case(seq_name, mers) =>
      mers foreach {case(word) =>
        if ( mer_counts contains word )
          mer_counts( word ) += 1
        else
          mer_counts += (word -> 1)
      }  // foreach
    }  // foreach
  }  // tallyMers
 
  // ***************************************************************************
  def bestMers() {
    val half = (proteins.keys.size + 1) / 2
    mer_counts foreach {case(word, count) =>
      if ( count >= half ) {
        top_mers += word
        // println( "count: " + count + ", half: " + half + ", word: " + word )
      }  // if
    }  // foreach

    // println( "bestMers: " + top_mers.mkString( "|" ) )
    // val positions = seq_tools.seqWords( proteins( name ).fasta.sequence, wordSize )
    // top_mers foreach {case(word) =>
    //   print( word + ":" + positions( word ).mkString( "," ) + " " )
    // }  // foreach
    // println
  }  // bestMers

  // ***************************************************************************
  def merCounts() {
    val seq = proteins( name ).fasta.sequence
    for ( i <- 0 until seq.length - (wordSize-1) ) {
      val kmer: String = seq.substring( i, i+wordSize )
      // if ( top_mers contains kmer )
        println( (i+1) + "\t" + kmer + "\t" + mer_counts( kmer ) )
    }  // for
    println
  }  // merCounts

  // ***************************************************************************
  def setAnchors() {
    val b = Map[String, ArrayBuffer[Int]]()		// y = mx + b for alignment
    val lookup = Map[String, Map[Int, Int]]()
    val half = (proteins.keys.size + 1) / 2
    val orders = Map[String, ArrayBuffer[Int]]()
    proteins foreach {case (align_name, protein) => 
      b += ( align_name -> ArrayBuffer[Int]() )
      orders += ( align_name -> ArrayBuffer[Int]() ) 
      lookup += ( align_name -> Map[Int, Int]() )
    }  // foreach

    // Traverse the words in the query.
    val seq = proteins( name ).fasta.sequence
    for ( i <- 0 until seq.length - (wordSize-1) ) {
      val word = seq.substring( i, i+wordSize )
      proteins foreach {case (align_name, protein) => 
        if ( name != align_name ) {
          if ( ( common_mers( align_name ) contains word ) && ( mer_counts( word ) >= half ) ) {
            if ( ( protein.positions contains word ) && ( protein.positions( word ).size == 1 ) ) {
              val word_pos = protein.positions( word ).last
              orders( align_name ) += word_pos
              lookup( align_name ) += (word_pos -> (i+1))
              if ( top_mers contains word )
                b( align_name ) += (word_pos - (i+1))
            }  // if
          }  // if
        }  // if
      }  // foreach
    }  // for

    // Initialize the alignment.
    align += name -> Map[Int, Char]()

    for ( i <- 0 until seq.length )
      align( name ) += (i+1) -> seq.charAt( i )
   
    orders foreach {case (align_name, order_list) =>
      val align_seq = proteins( align_name ).fasta.sequence

      if ( b( align_name ).size > 3 ) {
        val order = new Ordered( order_list )
        order.findBlocks()
        val linear = order.filter()

        val b_sorted = b( align_name ).sorted
        val b_target = b_sorted( b_sorted.size / 2 )
        // val pass0 = increasing( order_list, lookup(align_name) )
        val pass1 = bSelect( linear, lookup(align_name), b_target )
        // val pass3 = orderSelect( pass1 )

        val pass2 = increasing2( pass1, lookup(align_name) )

        pass2 foreach {case (i) => 
          val word2 = if ( i+wordSize-1 < align_seq.size ) align_seq.substring( i-1, i+wordSize-1 ) else "&#&" 
          val j = lookup( align_name )( i )
          // val m = i.toFloat / j.toFloat
          val b = j - i
        }  // foreach
        if ( pass2.size > 3 ) {
          align += align_name -> Map[Int, Char]()
          setWords( align_name, pass2, seq, lookup )
        }  // if
      }  // if
    }  // foreach

    // Determine the conservation levels.
    determineConservation()

    // Check the last position.
    snapAlign( seq )
  }  // setAnchors

  // ***************************************************************************
  def setWords( align_name: String, list: ArrayBuffer[Int], seq: String, lookup: Map[String, Map[Int, Int]] ) {
    if ( list.size < 3 )  return
 
    val align_seq = proteins( align_name ).fasta.sequence
 
    // Align sequence first segment.
    val align_start = list( 0 ) - 1
    val seq_start = lookup( align_name )( list( 0 ) ) - 1

    align5( align_name, align_seq, align_start, seq_start )

    // println( "setWords: " + align_name + " " + list.mkString( "|" ) )

    // Add the internal positions.
    for ( i <- 0 until list.size ) {
      val align_pos = list( i )
      val seq_pos = lookup( align_name )( align_pos )
      // println( "map: " + align_pos + " to " + seq_pos )
      for ( j <- 0 until 3 ) {
        align( align_name ) += (seq_pos+j) -> align_seq.charAt( align_pos+j-1 )
        // println( "--> " + align_name + " @ " + (align_pos+j) + " " + align_seq.charAt( align_pos+j-1 ) + " to " + (seq_pos+j) )
      }  // for
    }  // for

    // Align the gap positions.
    for ( i <- 0 until list.size-1 ) {
      val pos1 = list( i )
      val pos2 = list( i+1 )
      val seq_pos1 = lookup( align_name )( pos1 )
      val seq_pos2 = lookup( align_name )( pos2 )
      if ( ( pos2 > pos1+3 ) && ( seq_pos2 > seq_pos1+3 ) )
        subAlign( align_name, align_seq, pos1+3, pos2-1, seq, seq_pos1+3, seq_pos2-1 )
    }  // for

    // Align the sequence last segment.
    align3( align_name, align_seq, list.last+1, lookup( align_name )( list.last )+1 ) 
  }  // setWords

  // ***************************************************************************
  // Function checks if tuple is near neighbor tuples.
  def nearBy( align1: Int, align2: Int, align3: Int ): Boolean = {
    val dist1 = align2 - align1
    val dist2 = align3 - align2
    return ( (dist1 < 10) && (dist2 < 10) )
  }  // nearBy

  // ***************************************************************************
  def slope( x1: Int, y1: Int, x2: Int, y2: Int ): Float = {
    return (y2.toFloat - y1.toFloat) / (x2.toFloat - x1.toFloat)
  }  // slope

  // ***************************************************************************
  def goodSlope( m1: Float ): Boolean = {
    if ( ( 0.75F <= m1 ) && ( m1 <= 1.25F ) )
      return true
    return false
  }  // goodSlope

  // ***************************************************************************
  def similarSlopes( m1: Float, m2: Float ): Boolean = {
    if ( m1 * m2 < 0.0F )  return false

    return ( goodSlope( m1 ) && goodSlope( m2 ) )
  }  // similarSlopes

  // ***************************************************************************
  def bSelect( list: ArrayBuffer[Int], seqPos: Map[Int, Int], b_target: Int ): ArrayBuffer[Int] = {
    val inc = ArrayBuffer[Int]()

    if ( list.size < 3 )
      return inc

    for ( i <- 0 until list.size ) {
      val b = list( i ) - seqPos( list(i) )
      if ( ( b >= b_target - 15 ) && ( b <= b_target + 15 ) )
        inc += list( i )
    }  // for

    return inc
  }  // bSelect

  // ***************************************************************************
  def orderSelect( list: ArrayBuffer[Int] ): ArrayBuffer[Int] = {
    val inc = ArrayBuffer[Int]()

    if ( list.size < 3 )
      return inc

    for ( i <- 0 until list.size-2 ) {
      val last = if (i + 3 < list.size-1 ) i+3 else list.size-1
      var keep = true
      for ( j <- i+1 until last )
        if ( list( i ) > list( j ) )
          keep = false

      if ( keep )
        inc += list( i )
    }  // for

    // Add the last two elements.
    for ( i <- list.size-2 until list.size ) {
      if ( list( i ) > list( i-1 ) )
        inc += list( i )
    }  // for

    inc
  }  // orderSelect

  // ***************************************************************************
  def increasing( list: ArrayBuffer[Int], seqPos: Map[Int, Int] ): ArrayBuffer[Int] = {
    val inc = ArrayBuffer[Int]()

    if ( list.size < 3 )
      return inc

    // Add the first element of the list if ordered.
    // var near = nearBy( list(0), list(1), list(2) ) && nearBy( seqPos( list(0) ), seqPos( list(1) ), seqPos( list(2) ) )
    val m = slope( list(0), seqPos( list(0) ), list(1), seqPos( list(1) ) )
    // if ( ( list( 0 ) < list( 1 ) ) && ( list( 0 ) < list( 2 ) ) && goodSlope( m ) )
    if ( ( list( 0 ) < list( 1 ) ) && goodSlope( m ) )
      inc += list( 0 )

    // Scan the list.
    for ( i <- 1 until list.size-1 ) {
      // Don't keep words with large gap jumps.
      val m1 = slope( list(i-1), seqPos( list(i-1) ), list(i), seqPos( list(i) ) )
      val m2 = slope( list(i), seqPos( list(i) ), list( i+1 ), seqPos( list (i+1) ) )
      val good = similarSlopes( m1, m2 )
      // println( "slopes: " + list(i-1) + ":" + seqPos(list(i-1)) + "-" + list(i) + ":" + seqPos(list(i)) + "-" + list(i+1) + ":" + seqPos(list(i+1)) + " m1: " + m1 + " m2: " + m2 + " is good: " + good )
      if ( ( list( i-1) < list( i ) ) && ( list( i ) < list( i+1 ) ) && good )
        inc += list( i )
    }  // for

    // Check the last element.
    if ( ( inc.size > 0 ) && ( list.last > inc.last ) )
      inc += list.last

    return inc
  }  // increasing

  // ***************************************************************************
  def increasing2( list: ArrayBuffer[Int], seqPos: Map[Int, Int] ): ArrayBuffer[Int] = {
    val inc = ArrayBuffer[Int]()

    if ( list.size < 3 )
      return inc

    // Add the first element of the list if ordered.
    val m = slope( list(0), seqPos( list(0) ), list(1), seqPos( list(1) ) )
    // if ( ( list( 0 ) < list( 1 ) ) && ( list( 0 ) < list( 2 ) ) && goodSlope( m ) )
    if ( ( list( 0 ) < list( 1 ) ) && goodSlope( m ) )
      inc += list( 0 )

    // Scan the list.
    for ( i <- 1 until list.size-1 ) {
      // Don't keep words with large gap jumps.
      val m1 = slope( list(i-1), seqPos( list(i-1) ), list(i), seqPos( list(i) ) )
      val m2 = slope( list(i), seqPos( list(i) ), list( i+1 ), seqPos( list (i+1) ) )
      val good = similarSlopes( m1, m2 )
      // println( "slopes: " + list(i-1) + ":" + seqPos(list(i-1)) + "-" + list(i) + ":" + seqPos(list(i)) + "-" + list(i+1) + ":" + seqPos(list(i+1)) + " m1: " + m1 + " m2: " + m2 + " is good: " + good )
      if ( ( list( i-1) < list( i ) ) && ( list( i ) < list( i+1 ) ) && good )
        inc += list( i )
    }  // for

    // Check the last element.
    if ( ( inc.size > 0 ) && ( list.last > inc.last ) )
      inc += list.last

    return inc
  }  // increasing2

  // ***************************************************************************
  def findGap( align_seq: String, align1: Int, align2: Int, seq: String, seq1: Int, seq2: Int ): Tuple2[Int, Int] = {
    val gap_size = (seq2 - seq1) - (align2 - align1)

    var best_start = seq1
    var best_count = 0
    for ( gap_start <- seq1 until (seq2-gap_size+1) ) {
      var count = 0
      for ( i <- align1 until align2+1 ) {
        if ( i-align1+seq1 < gap_start ) {
          if ( align_seq.charAt( i-1 ) == seq.charAt( i-align1+seq1-1 ) ) 
            count += 1
        }   
        else {
          if ( align_seq.charAt( i-1 ) == seq.charAt( i-align1+seq1+gap_size-1 ) ) 
            count += 1
        }  // if
      }  // for

       if ( count > best_count ) {
         best_count = count
         best_start = gap_start
       }  // if
    }  // for

    return (best_start, best_start+gap_size-1)
  }  // findGap

  // ***************************************************************************
  def findInsert( align_seq: String, align1: Int, align2: Int, seq: String, seq1: Int, seq2: Int ): Tuple2[Int, Int] = {
    val ins_size = (align2 - align1) - (seq2 - seq1)

    var best_start = align1
    var best_count = 0
    for ( ins_start <- align1 until (align2-ins_size+1) ) {
      var count = 0
      for ( i <- seq1 until seq2+1 ) {
        if ( i-seq1+align1 < ins_start ) {
          if ( align_seq.charAt( i-seq1+align1-1 ) == seq.charAt( i-1 ) ) 
            count += 1
        }   
        else {
          if ( align_seq.charAt( i-seq1+align1+ins_size-1 ) == seq.charAt( i-1 ) ) 
            count += 1
        }  // if
      }  // for

       if ( count > best_count ) {
         best_count = count
         best_start = ins_start
       }  // if
    }  // for

    return (best_start, best_start+ins_size-1)
  }  // findInsert

  // ***************************************************************************
  def subAlign( align_name: String, align_seq: String, align1: Int, align2: Int, seq: String, seq1: Int, seq2: Int ) {
    // println( "subAlign: " + align_name + " [" + align1 + "-" + align2 + "] X [" + seq1 + "-" + seq2 + "]" )

    // Check for matching lengths and a segment.
    val sub_sequence = align_seq.substring( align1-1, align2 )
    if ( align2 - align1 == seq2 - seq1 ) {
      for ( i <- 0 until (seq2-seq1+1) ) {
        val aa = if ( align1+i < align_seq.size ) align_seq.charAt( align1+i-1 ) else '-'
        align( align_name ) += ( seq1+i -> aa )
      }  // for

      if ( ( segments contains align_name ) == false )
        segments += align_name -> Map[Int, Insert]()
      segments( align_name ) += (seq1) -> new Insert( align_name, seq1-1, align_seq.substring( align1-1, align2 ), new Window( align1, align2, sub_sequence ) )
    }  // if
    else
      if ( align2 - align1 > seq2 - seq1 ) {
        // Optimize subalignment and add insert of extra residues.
        // println( "subAlign: addInsert " + align_name + " [" + align1 + "-" + align2 + "] X [" + seq1 + "-" + seq2 + "]" )

        // ins_start and ins_end are based on the align_seq coordindates.
        val (ins_start, ins_end) = findInsert( align_seq, align1, align2, seq, seq1, seq2 )
        // println( "###Insert: " + align_name + " [" + align1 + ":" + align2 + "] vs. [" + seq1 + ":" + seq2 + "] ins [" + ins_start + "-" + ins_end + "]" )
        if ( ( inserts contains align_name ) == false )
          inserts += align_name -> Map[Int, Insert]()
        inserts( align_name ) += (ins_start-1) -> new Insert( align_name, seq1+ins_start-align1-1, align_seq.substring( ins_start-1, ins_end ), new Window( align1, align2, sub_sequence ) )

        // Align before the insert.
        if ( align1 < ins_start ) 
          for ( i <- align1 until ins_start )
            align( align_name ) += i-align1+seq1 -> align_seq.toLowerCase().charAt( i-1 )

        // Align after the insert.
        if ( ins_end < align2 )
          for ( i <- ins_end+1 until align2+1 )
            align( align_name ) += seq2-(align2-i) -> align_seq.toLowerCase().charAt( i-1 )
      }
      else {
        // Optimize subalignment and add gap.
        // gap_start and gap_end are based on the seq alignment coordinates.
        val (gap_start, gap_end) = findGap( align_seq, align1, align2, seq, seq1, seq2 )
        // println( "###Gap: " + align_name + " [" + align1 + ":" + align2 + "] vs. [" + seq1 + ":" + seq2 + "] gap [" + gap_start + "-" + gap_end + "]" )
        if ( ( gaps contains align_name ) == false )
          gaps += align_name -> Map[Int, Gap]()
        gaps( align_name ) += gap_start -> new Gap( align_name, gap_start, gap_end, new Window( align1, align2, sub_sequence ) )

        // Align before the gap.
        if ( gap_start > seq1 ) {
          for ( i <- seq1 until gap_start)
            align( align_name ) += i -> align_seq.toLowerCase().charAt( i-seq1+align1-1 )
        }  // if

        // Set the gap characters.
        for ( i <- gap_start until gap_end+1 )
          align( align_name ) += i -> '-'

        // Align after the gap. 
        if ( gap_end < seq2 ) {
          for ( i <- gap_end+1 until seq2+1 )
            align( align_name ) += i -> align_seq.toLowerCase().charAt( align2-(seq2-i)-1 )
        }  // if
      }  // if
  }  // subAlign

  // ***************************************************************************
  def align3( align_name: String, align_seq: String, align_end: Int, seq_end: Int ) {
    if ( align_end >= align_seq.size )  return

    for ( i <- 0 until align_seq.size-align_end ) {
      align( align_name ) += ( seq_end+i+1 -> align_seq.toLowerCase().charAt( align_end+i ) )
    }  // for
  }  // align3

  // ***************************************************************************
  def align5( align_name: String, align_seq: String, align_start: Int, seq_start: Int ) {
    // println( "align5: " + align_name + " seq_start: " + seq_start + ", align_start: " + align_start )

    if ( align_start < 1 ) return

    // Case 1: same length
    if ( align_start == seq_start ) {
      for ( i <- 0 until align_start )
        align( align_name ) += ( i+1 -> align_seq.charAt( i ) )
      return
    }  // if
   
    // Case 2: align sequence is shorter.
    if ( align_start < seq_start ) {
      for ( i <- 0 until align_start )
        align( align_name ) += ( seq_start-align_start+i+1 -> align_seq.charAt( i ) )
      return
    }  // if

    // Case 3: align sequence is longer.
    for ( i <- 0 until seq_start )
      align( align_name ) += ( i+1 -> align_seq.charAt( align_start-seq_start+i ) )

    if ( ( inserts contains align_name ) == false )
      inserts += align_name -> Map[Int, Insert]()

    inserts( align_name ) += 0 -> new Insert( align_name, 0, align_seq.slice( 0, align_start-seq_start ), new Window( 0, 0, "" ) )
  }  // align5

  // ***************************************************************************
  def groupLoops() {
    if ( inserts.size < 1 )  return

    // Add all of the insertions to loops. 
    inserts foreach { case(align_name, inser) =>
      inser foreach { case(align1, ins) => 
        if ( ( ins.after > 0 ) && ( ins.mark == false ) ) {
          val loop = new Loop( ins.after-5, ins.after+5, new ArrayBuffer[Gap](), new ArrayBuffer[Insert]() )
          ins.mark = true
          loop.inserts += ins
          loop.names += align_name
  
          // Check for overlapping inserts
          inserts foreach { case(aln_name, insertList) =>
            insertList foreach { case(align1, insert) => 
              if ( ( insert.after > 0 ) && ( insert.mark == false ) && ( insert.after < ins.after+5 ) && ( insert.after > ins.after-5 ) ) {
                loop.inserts += insert
                insert.mark = true
                loop.names += aln_name
              }  // if
            }  // foreach
          }  // foreach
  
          // Check for overlapping gaps.
          gaps foreach { case(gapName, gapList) =>
            gapList foreach { case(pos, gap) =>
              if ( ( gap.mark == false ) && ( gap.start < ins.after+5 ) && ( gap.end > ins.after-5 ) ) {
                gap.mark = true
                loop.gaps += gap
                loop.names += gapName
              }  // if
            }  // foreach
          }  // foreach

          loop.findSite()
          loops += loop.after -> loop
        }  // if    
      }  // foreach
    }  // foreach

    // Create loops for the remaining unmarked gaps.
    gaps foreach { case(gapName, gapList) =>
      gapList foreach { case(pos, gap) =>
        if ( gap.mark == false ) {
          // println( "*** unmarked gap: " + gapName + " " + gap.start + "-" + gap.end )
          val loop = new Loop( gap.start-5, gap.end+5, new ArrayBuffer[Gap](), new ArrayBuffer[Insert]() )
          gap.mark = true
          loop.gaps += gap
          loop.names += gapName

          // Check for overlapping gaps.
          gaps foreach { case(gap_name, gap_list) =>
            gap_list foreach { case(pos2, gap2) =>
              if ( ( gap2.mark == false ) && ( gap2.start < gap.end+5 ) && ( gap2.end > gap.start-5 ) ) {
                gap2.mark = true
                loop.gaps += gap2
                loop.names += gap_name
              }  // if
            }  // foreach
          }  // foreach         

          loop.findGapSite()
          loops += loop.after -> loop
        }  // if
      }  // foreach
    }  // foreach
  }  // groupLoops

  // ***************************************************************************
  def largestPrefix(): Int = {
    var size = 0
    inserts foreach {case(alignName, ins) =>
      if ( ( ins contains 0 ) && ( ins( 0 ).seq.size > size ) )
        size = ins( 0 ).seq.size
    }  // foreach

    return size
  }  // largestPrefix

  // ***************************************************************************
  def rightAdjust( preSeq: String, preSize: Int ): String = {
    val str = new StringBuilder()
    val dotCount = preSize - preSeq.size
    if ( dotCount > 0 )
      for ( i <- 0 until dotCount )  str += '.'

    str.mkString + preSeq
  }  // rightAdjust

  // ***************************************************************************
  def findCore( win_start: Int, win_end: Int, align_name: String, core_seq: String ): Tuple2[Int, Int] = {
    var best_count = 0
    var best_start = -1
    for ( i <- win_start until win_end-core_seq.size ) {
      var count = 0
      for ( j <- 0 until core_seq.size ) {
        if ( ( align( align_name) contains i+j ) && ( core_seq.charAt( j ) == align( align_name )( i+j ).toUpper ) ) {
          count += 1
        }  // if
      }  // for

      if ( count > best_count ) {
        best_start = i
        best_count = count
      }  // if
    }  // for

    new Tuple2(best_start, best_count)
  }  // findCore

  // ***************************************************************************
  def sliceAlign( align_name: String, start: Int, end: Int ): String = {
    val str = new StringBuilder()
    for ( i <- start to end )
      str += align( align_name )( i )
    str.toString
  }  // sliceAlign

  // ***************************************************************************
  def adjustAlignment( align_name: String, win_start: Int, win_end: Int, block_start: Int, block_end: Int, core_start: Int, core_seq: String ) {
    if ( core_start < block_start ) {
      // Slide alignment to the right
      val delta = block_start - core_start 

      // Set up the new insert.
      val ins_seq = sliceAlign( align_name, block_end - delta, block_end )
      if ( ( inserts contains align_name ) == false )
        inserts += align_name -> Map[Int, Insert]()
      inserts( align_name ) += block_end -> new Insert( align_name, block_end, ins_seq, new Window( block_end+1, win_end, ins_seq ) )

      // Shift the alignment.
      for ( i <- block_end to block_start by -1 ) 
        align( align_name )( i ) = align( align_name )( i-delta )

      // Back fill with gaps.
      for ( i <- core_start until core_start+delta )
        align( align_name )( i ) = '-'

      // Create a new gap data structure.
      if ( ( gaps contains align_name ) == false )
        gaps += align_name -> Map[Int, Gap]()
      gaps( align_name ) += core_start -> new Gap( align_name, core_start, core_start+delta, new Window( win_start, core_start+delta, "" ) )
    } else {
      // Slide alignment to the left
      val delta = core_start - block_start 
      val ins_seq = sliceAlign( align_name, block_start, block_start + delta )
      if ( ( inserts contains align_name ) == false )
        inserts += align_name -> Map[Int, Insert]()
      inserts( align_name ) += (block_start-1) -> new Insert( align_name, block_start-1, ins_seq, new Window( win_start, block_start-1, ins_seq ) )

      // Shift the alignment.
      for ( i <- block_start to block_end )
        align( align_name )( i ) = align( align_name)( i+delta )

      // Fill in the new gaps.
      for ( i <- block_end+1 to block_end+delta )
        align( align_name )( i ) = '-'

      // Create a new gap data structure.
      if ( ( gaps contains align_name ) == false )
        gaps += align_name -> Map[Int, Gap]()
      gaps( align_name ) += (block_end+1) -> new Gap( align_name, block_end+1, block_end+delta, new Window( block_end+1, win_end, "" ) )
    }  // if
  }  // adjustAlignment

  // ***************************************************************************
  def optimizeBlock( win_start: Int, win_end: Int, block_start: Int, block_end: Int, core_seq: String ) {
    // println( "optimizeBlock " + win_start + "-" + win_end + " " + core_seq )
    align foreach {case(align_name, seq_align) =>
      val (core_start, core_count) = findCore( win_start, win_end, align_name, core_seq )
      val core_seg = new StringBuilder()
      for ( i <- 0 until core_seq.size )
        if ( align( align_name ) contains core_start+i )
          core_seg += align( align_name )( core_start+i )
        else '?'
      if ( ( core_start != block_start ) && ( core_start > 0 ) && ( ( core_count >= 4 ) || ( ( core_count == 3 ) && ( core_seq.size == 5 ) ) ) ) {
        // println( "  " + align_name + "  core: " + core_start + " block: " + block_start + "  " + core_seg + "  " + core_count )
        adjustAlignment( align_name, win_start, win_end, block_start, block_end, core_start, core_seq )
      }  // if
    }  // foreach

  }  // optimzeBlock

  // ***************************************************************************
  def optimizeCore( core_blocks: ArrayBuffer[Window] ) {
    if ( core_blocks.size > 0 ) {
      var win_start = if ( core_blocks( 0 ).start >= 10 ) core_blocks( 0 ).start-10 else 1
      for ( i <- 0 until core_blocks.size ) {
        var win_end = if ( i < core_blocks.size-1 ) core_blocks( i+1 ).start-1 else core_blocks( i ).end + 10
        win_end = if ( win_end < align( name ).size ) win_end else align( name ).size
        optimizeBlock( win_start, win_end, core_blocks( i ).start, core_blocks( i ).end, core_blocks( i ).seq )
        win_start = core_blocks( i ).end + 1
      }  // for
    }  // if
  }  // optimizeCore 

  // ***************************************************************************
  def copyBlock( full_end: Int, align_end: Int, block_end: Int ) {
    // println( "*** copyBlock: from " + (align_end+1) + " to " + block_end )

    align foreach {case(align_name, seq_align) =>
      for ( i <- align_end+1 to block_end )
        if ( align( align_name ) contains i )
          full_align( align_name )(full_end+i-align_end ) = align( align_name )( i )
        else
          full_align( align_name )(full_end+i-align_end ) = '.'
    }  // foreach

    // snapFullAlign()
  }  // copyBlock

  // ***************************************************************************
  def optimizeMsa() {
    // println
    // println( "Optimze MSA:" )

    // Set up the Pre-alignment sequences.
    val pre_size = largestPrefix()

    align foreach {case(align_name, seq_align) =>
      val pre_sequence = if ( ( inserts contains align_name) && ( inserts( align_name ) contains 0 ) ) rightAdjust( inserts( align_name )( 0 ).seq, pre_size ) else rightAdjust( "", pre_size )

      full_align += align_name -> Map[Int, Char]()
      if ( pre_size > 0 )
        for ( i <- 0 until pre_size )  full_align( align_name ) += i+1 -> pre_sequence.charAt( i )
    }  // foreach

    var full_end = pre_size
    var align_end = 0

    val loop_starts = loops.keys.toList.sorted.foreach{case(start) =>
      val loop = loops( start )
      copyBlock( full_end, align_end, loop.after )

      // println( "full_end.0: " + full_end )
      full_end += loop.after - align_end
      // println( "full_end.1: " + full_end + ", loop.after: " + loop.after + ", align_end: " + align_end + ", loop.size: " + loop.size )
      align_end = loop.after

      // println( "loop: " + loop.start + "-" + loop.end + " after: " + loop.after + ", size: " + loop.size )
      loop.inserts foreach { case(insert) => 
        // println( "  " + insert.name+ " insert: " + insert.after + " " + insert.seq + " in [" + insert.window.start + "-" + insert.window.end + "]  " + insert.window.seq ) 
        for ( i <- 0 until insert.seq.size )
          full_align( insert.name )( full_end+i+1 ) = insert.seq.charAt( i )
        if ( insert.seq.size < loop.size )
        for ( i <- insert.seq.size+1 to loop.size )
          full_align( insert.name )( full_end+i ) = '-'
      }  // foreach

      loop.gaps foreach {case (gap) => 
          val gap_size = gap.end - gap.start + 1
          // println( "  " + gap.name + " gap: (" + gap_size + ") " + gap.start + "-" + gap.end + " in [" + gap.window.start + "-" + gap.window.end + "]  " + gap.window.seq ) 
          for ( i <- 1 to loop.size )
            full_align( gap.name )( full_end+i ) = '-'
      }  // foreach

      // Add alignment gaps to sequences not included in this loop.
      if ( loop.size > 0 )
        align foreach {case(aln_name, sq_align) =>
          if ( loop.names.contains( aln_name ) == false )
            for ( i <- 1 to loop.size )
              full_align( aln_name )( full_end+i ) = '-'
        }  // foreach

      full_end += loop.size
    }  // foreach

    // Copy the MSA after the last loop.
    var max_key = 0
    align foreach {case(align_name, seq_align) => 
      val last_key = align( align_name ).keys.toList.sorted.last
      if ( last_key > max_key )  max_key = last_key
    }  // foreach

    copyBlock( full_end, align_end, max_key )

    // snapFullAlign()
  }  // optimizeMsa

  // ***************************************************************************
  def snapAlign( seq: String ) = {
    println( "MSA:" ) 
    val preSize = largestPrefix()

    groupLoops()

    gene_names foreach {case(gene_name, count) =>
      align foreach {case(align_name, seq_align) =>
        val gene = if ( proteins( align_name ).fasta.annotation contains "gene" ) proteins( align_name ).fasta.annotation( "gene" ).toUpperCase() else "other"

        if ( gene == gene_name ) {
          print( align_name + "\t" )
          // if ( ( inserts contains align_name) && ( inserts( align_name ) contains 0 ) )
          //   print( rightAdjust( inserts( align_name )( 0 ).seq, preSize ) )
          // else
          //   print( rightAdjust( "", preSize ) )
    
          if ( preSize > 0 )  print( "|" )
    
          for ( i <- 1 until seq.size+1 )
            if ( align( align_name ) contains i )
              print( align( align_name )( i ) )
            else
              print( "." )
          println( "  " + align_name + "  " + gene_name )
        }  // if
      }  // foreach
    }  // foreach
  } // snapAlign

  // ***************************************************************************
  def snapFullAlign() = {
    println( "Full MSA:" ) 

    gene_names foreach {case(gene_name, count) =>
      full_align foreach {case(align_name, seq_align) =>
        val gene = if ( proteins( align_name ).fasta.annotation contains "gene" ) proteins( align_name ).fasta.annotation( "gene" ).toUpperCase() else "other"

        if ( gene == gene_name ) {
          print( align_name + "\t" )
    
          val align_end = full_align( align_name ).keys.toList.sorted.last 

          // println( "snapFullAlign: " + align_name + " last: " + align_end )

          for ( i <- 1 until align_end )
            if ( full_align( align_name ) contains i )
              print( full_align( align_name )( i ) )
            else
              print( "!" )
          println( "  " + align_name + "  " + gene_name )
        }  // if
      }  // foreach
    }  // foreach
  } // snapFullAlign
  
  // ***************************************************************************
  def sliceClique( readers: Map[String, PdbReader], minConservation: Int ) {
    // Slice the clique from each PDB file.
    readers foreach {case(pdb_name, reader) =>
      if ( align contains pdb_name ) {
        // println( "sliceClique: pdb " + pdb_name )
        val tokens = pdb_name.split( "_" )
        val chain: Char = tokens( 1 ).charAt( 0 )
        val molecule = reader.residues( chain )
        val out = new OutputFile( pdb_name + ".clique" )
    
        var j: Int = 0
        for ( i <- 1 until (align(name).size+1) ) {
          if ( ( conservation contains i ) && ( conservation( i ) >= minConservation ) ) 
            // molecule( j ).atoms foreach {case(atom) => out.write( "j: " + j + ", i: " + i + "  " + align( pdb_name )( i ) + " pdb: " + atom.letter + " " + atom.line + "\n" ) }
            molecule( j ).atoms foreach {case(atom) => out.write( atom.line + "\n" ) }
          // else
            // if ( ( conservation contains i ) && ( align( pdb_name ) contains i ) && ( conservation( i ) < minConservation ) && ( conservation( i ) > 10 ) ) 
              // molecule( j ).atoms foreach {case(atom) => out.write( "--- j: " + j + ", i: " + i + "  " + align( pdb_name )( i ) + " pdb: " + atom.letter + " " + atom.line + "\n" ) }
    
          if ( ( align( pdb_name ) contains i ) && ( align( pdb_name )( i ) != '.' ) )
            j += 1 
        }  // for
      }  // if
    }  // foreach
  }  // sliceClique
 
  // ***************************************************************************
   def snapLoops {
     println
     println( "Loops:" )
     loops foreach { case(start, loop) =>
       println( "loop: " + loop.start + "-" + loop.end + " after: " + loop.after + ", size: " + loop.size )
       loop.inserts foreach { case(insert) => println( "  " + insert.name+ " insert: " + insert.after + " " + insert.seq + " in [" + insert.window.start + "-" + insert.window.end + "] " ) }
       loop.gaps foreach {case (gap) => println( "  " + gap.name + " gap:" + gap.end + " in [" + gap.window.start + "-" + gap.window.end + "] " ) }
     }  // foreach


   }  // snapLoops 

  // ***************************************************************************
  def getAlignment( align_name: String ): String = {
    val align = new StringBuilder()
    val align_end = full_align( align_name ).keys.toList.sorted.last 
    for ( i <- 1 until align_end )
      if ( full_align( align_name ) contains i )
        align += full_align( align_name )( i )

    align.toString
  }  // getAlignment 

  // ***************************************************************************
  def msaName : String = {
    var best_count = 0
    var best_name = "other"
    gene_names foreach {case(gene_name, count) =>
      if ( count > best_count ) {
        best_name = gene_name
        best_count = count
      }  // if
    }  // foreach

    best_name.replaceAll(" ", "").filterNot( _ == '/' )
  }  // msaName

  // ***************************************************************************
  def writeMsa( name: String ) {
    val msa_name = msaName
    System.out.println( "writeMsa: " + msa_name + ", name: " + name )
    val msa_file = new OutputFile( msa_name + ".fa" )
    full_align foreach {case(align_name, seq_align) =>
      val gene = proteins( align_name )
      msa_file.write( ">" + gene.fasta.name + " " + gene.fasta.description + "\n" + gene.fasta.toBlock( getAlignment( align_name ) ) )
    }  // foreach
  }  // writeMsa

  // ***************************************************************************
}  // class Msa
