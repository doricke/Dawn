
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

// *****************************************************************************
class SubAlignment( val target_seq: String, val align_name: String, val align_seq: String, wordSize: Int ) {

  val seq_tools = new SeqTools()

  val align_mers: collection.mutable.Set[String] = seq_tools.seqSplit( align_seq, wordSize )

  val align_pos: collection.mutable.Map[String, ArrayBuffer[Int]] = seq_tools.seqWords( align_seq, wordSize )

  val gaps = Map[String, Map[Int, Gap]]()		// gap start-end

  val inserts = Map[String, Map[Int, Insert]]() 	// after:insertion

  val align_segment = Map[Int, Char]()			// alignment segment

  val target_mers: collection.mutable.Set[String] = seq_tools.seqSplit( target_seq, wordSize )

  val target_pos: collection.mutable.Map[String, ArrayBuffer[Int]] = seq_tools.seqWords( target_seq, wordSize )

  // println( "SubAlignment: " + align_name + "  " + target_seq + " X " + align_seq )

  // ***************************************************************************
  def setAnchors() {
    val b = ArrayBuffer[Int]()		// y = mx + b for alignment
    val lookup = Map[Int, Int]()
    val orders = ArrayBuffer[Int]()

    // Traverse the words in the query.
    for ( i <- 0 until target_seq.length - (wordSize-1) ) {
      val word = target_seq.substring( i, i+wordSize )
      // println( "word: " + word + " " + i )
      if ( align_mers contains word )
        if ( ( align_pos contains word ) && ( align_pos( word ).size == 1 ) ) {
          val word_pos = align_pos( word ).last
          orders += word_pos
          lookup += (word_pos -> (i+1))
          b += (word_pos - (i+1))
          // println( "lookup: " + align_name + " " + word + " " + word_pos + " @ " + (i+1) + " " )
        }  // if
    }  // for
   
    // println
    // println( align_name + "  pre-pass: " + order_list.mkString( "|" ) )
    if ( b.size > 3 ) {
      val order = new Ordered( orders )
      order.findBlocks()
      val linear = order.filter()
      println( align_name + "  filtered: " + linear.mkString( "|" ) )

      val b_sorted = b.sorted
      val b_target = b_sorted( b_sorted.size / 2 )
      // val pass0 = increasing( order_list, lookup )
      // println( align_name + "  pass0: " + pass0.mkString( "|" ) )
      val pass1 = bSelect( linear, lookup, b_target )
      println( align_name + "  pass1: " + pass1.mkString( "|" ) )
      // val pass3 = orderSelect( pass1 )
      // println( align_name + "  pass1b: " + pass3.mkString( "|" ) )

      val pass2 = increasing2( pass1, lookup )
      println( align_name + "  pass2: " + pass2.mkString( "|" ) )

      // print( "pass2: " + align_name + " " )
      pass2 foreach {case (i) => 
        val word2 = if ( i+wordSize-1 < align_seq.size ) align_seq.substring( i-1, i+wordSize-1 ) else "&#&" 
        val j = lookup( i )
        // val m = i.toFloat / j.toFloat
        val b = j - i
        // print( word2 + ":" + i + "/" + j + "|" + b + " " ) 
      }  // foreach
      // println
      if ( pass2.size > 3 ) {
        // println( "@@@ setWords for " + align_name + ", words: " + pass2.size )
        setWords( pass2, lookup )
      }  // if
    }  // if
  }  // setAnchors

  // ***************************************************************************
  def setWords( list: ArrayBuffer[Int], lookup: Map[Int, Int] ) {
    if ( list.size < 3 )  return
 
    // Align sequence first segment.
    val align_start = list( 0 ) - 1
    val seq_start = lookup( list( 0 ) ) - 1

    // align5( align_name, align_seq, align_start, seq_start )

    // println( "setWords: " + align_name + " " + list.mkString( "|" ) )

    // Add the internal positions.
    for ( i <- 0 until list.size ) {
      val align_pos = list( i )
      val seq_pos = lookup( align_pos )
      // println( "map: " + align_pos + " to " + seq_pos )
      for ( j <- 0 until 3 ) {
        align_segment += (seq_pos+j) -> align_seq.charAt( align_pos+j-1 )
        // println( "--> " + align_name + " @ " + (align_pos+j) + " " + align_seq.charAt( align_pos+j-1 ) + " to " + (seq_pos+j) )
      }  // for
    }  // for

    // Align the gap positions.
    for ( i <- 0 until list.size-1 ) {
      val pos1 = list( i )
      val pos2 = list( i+1 )
      val seq_pos1 = lookup( pos1 )
      val seq_pos2 = lookup( pos2 )
      if ( ( pos2 > pos1+3 ) || ( seq_pos2 > seq_pos1+3 ) )
        subAlign( pos1+3, pos2-1, seq_pos1+3, seq_pos2-1 )
    }  // for

    // Align the sequence last segment.
    // align3( list.last+1, lookup( list.last )+1 ) 
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
  def findGap( align1: Int, align2: Int, seq1: Int, seq2: Int ): Tuple2[Int, Int] = {
    val gap_size = (seq2 - seq1) - (align2 - align1)

    var best_start = seq1
    var best_count = 0
    for ( gap_start <- seq1 until (seq2-gap_size+1) ) {
      var count = 0
      for ( i <- align1 until align2+1 ) {
        if ( i-align1+seq1 < gap_start ) {
          if ( align_seq.charAt( i-1 ) == target_seq.charAt( i-align1+seq1-1 ) ) 
            count += 1
        }   
        else {
          if ( align_seq.charAt( i-1 ) == target_seq.charAt( i-align1+seq1+gap_size-1 ) ) 
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
  def findInsert( align1: Int, align2: Int, seq1: Int, seq2: Int ): Tuple2[Int, Int] = {
    val ins_size = (align2 - align1) - (seq2 - seq1)

    var best_start = align1
    var best_count = 0
    for ( ins_start <- align1 until (align2-ins_size+1) ) {
      var count = 0
      for ( i <- seq1 until seq2+1 ) {
        if ( i-seq1+align1 < ins_start ) {
          if ( align_seq.charAt( i-seq1+align1-1 ) == target_seq.charAt( i-1 ) ) 
            count += 1
        }   
        else {
          if ( align_seq.charAt( i-seq1+align1+ins_size-1 ) == target_seq.charAt( i-1 ) ) 
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
  def subAlign( align1: Int, align2: Int, seq1: Int, seq2: Int ) {
    // println( "subAlign: " + align_name + " [" + align1 + "-" + align2 + "] X [" + seq1 + "-" + seq2 + "]" )

    // Check for matching lengths and a segment.
    val sub_sequence = align_seq.substring( align1-1, align2 )
    if ( align2 - align1 == seq2 - seq1 ) {
      for ( i <- 0 until (-seq1+1) ) {
        val aa = if ( align1+i < align_seq.size ) align_seq.charAt( align1+i-1 ) else '@'
        align_segment += ( seq1+i -> aa )
      }  // for
    }  // if
    else
      if ( align2 - align1 > seq2 - seq1 ) {
        // Optimize subalignment and add insert of extra residues.
        // println( "subAlign: addInsert " + align_name + " [" + align1 + "-" + align2 + "] X [" + seq1 + "-" + seq2 + "]" )

        // ins_start and ins_end are based on the align_seq coordindates.
        val (ins_start, ins_end) = findInsert( align1, align2, seq1, seq2 )
        // println( "###Insert: " + align_name + " [" + align1 + ":" + align2 + "] vs. [" + seq1 + ":" + seq2 + "] ins [" + ins_start + "-" + ins_end + "]" )
        if ( ( inserts contains align_name ) == false )
          inserts += align_name -> Map[Int, Insert]()
        inserts( align_name ) += (ins_start-1) -> new Insert( align_name, ins_start-1, align_seq.substring( ins_start-1, ins_end ), new Window( align1, align2, sub_sequence ) )

        // Align before the insert.
        if ( align1 < ins_start ) 
          for ( i <- align1 until ins_start )
            align_segment += i-align1+seq1 -> align_seq.toLowerCase().charAt( i-1 )

        // Align after the insert.
        if ( ins_end < align2 )
          for ( i <- ins_end+1 until align2+1 )
            align_segment += seq2-(align2-i) -> align_seq.toLowerCase().charAt( i-1 )
      }
      else {
        // Optimize subalignment and add gap.
        // gap_start and gap_end are based on the seq alignment coordinates.
        val (gap_start, gap_end) = findGap( align1, align2, seq1, seq2 )
        // println( "###Gap: " + align_name + " [" + align1 + ":" + align2 + "] vs. [" + seq1 + ":" + seq2 + "] gap [" + gap_start + "-" + gap_end + "]" )
        if ( ( gaps contains align_name ) == false )
          gaps += align_name -> Map[Int, Gap]()
        gaps( align_name ) += gap_start -> new Gap( align_name, gap_start, gap_end, new Window( align1, align2, sub_sequence ) )

        // Align before the gap.
        if ( gap_start > seq1 ) {
          for ( i <- seq1 until gap_start)
            align_segment += i -> align_seq.toLowerCase().charAt( i-seq1+align1-1 )
        }  // if

        // Set the gap characters.
        for ( i <- gap_start until gap_end+1 )
          align_segment += i -> '-'

        // Align after the gap. 
        if ( gap_end < seq2 ) {
          for ( i <- gap_end+1 until seq2+1 )
            align_segment += i -> align_seq.toLowerCase().charAt( align2-(seq2-i)-1 )
        }  // if
      }  // if
  }  // subAlign

  // ***************************************************************************
}  // class SubAlignment
