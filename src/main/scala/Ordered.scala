
import scala.collection.mutable.ArrayBuffer

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

// ****************************************************************************
class Ordered( val list: ArrayBuffer[Int] ) {

  val blocks: ArrayBuffer[Block] = new ArrayBuffer[Block]()

  val linear: ArrayBuffer[Int] = new ArrayBuffer[Int]()

// ****************************************************************************
  def findBlocks() {
    if ( list.size <= 3 )  return 

    var in_block = false
    var block_start = -1
    var block_end = -1

    for ( i <- 0 until list.size-1 ) {
      if ( in_block ) {
        if ( list(i)+1 == list(i+1) ) 
          block_end = i+1
        else {
          if ( block_end - block_start + 1 >= 3 )
            blocks += new Block( block_start, block_end, list( block_start ), list( block_end ) )

          block_start = -1
          block_end = -1
          in_block = false
        }  // if
      }
      else {
        if ( list(i)+1 == list(i+1) ) {
          in_block = true
          block_start = i
          block_end = i+1
        }  // if
      }  // if
    }  // for

    // println
    // println( "Blocks" )
    // blocks foreach { case(block) =>
    //   println( "block: " + block.start + "-" + block.end + "  [" + block.list_start + "-" + block.list_end + "]" )
    // }  // foreach
  }  // findBlocks

// ****************************************************************************
  def filter(): ArrayBuffer[Int] = {

    // Traverse the blocks filtering the list.
    var start = 0
    for ( i <- 0 until blocks.size ) {
      // Filter the list before this block.
      if ( start < blocks( i ).start ) {
        for ( j <- start until blocks( i ).start )
          if ( list( j ) < blocks( i ).list_start )
            linear += list( j )
      }  // if

      // Add this block
      for ( j <- blocks( i ).list_start until blocks( i ).list_end )
        linear += j

      start = blocks( i ).end + 1
    }  // for

    // Add after the last block.
    if ( blocks.size > 0 ) {
      val last_block = blocks( blocks.size-1 )
      if ( last_block.end+1 < list.size ) {
        for ( j <- last_block.end+1 until list.size )
          if ( list( j ) > last_block.list_end )
            linear += list( j )
      }  // if
    }  // if

    return linear 
  }  // filter

  // **************************************************************************
}  // class Ordered
