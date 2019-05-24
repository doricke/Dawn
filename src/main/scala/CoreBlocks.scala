
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
class CoreBlocks( val consensus: String ) {

  val blocks = new ArrayBuffer[Window]()

// *****************************************************************************
  def makeBlocks() {
    var block_start = -1
    var block_end = -1
    var in_block = false
    for ( i <- 0 until consensus.size ) {
      if ( consensus.charAt( i ) != ' ' ) {
        if ( in_block == false ) {
          in_block = true
          block_start = i
        }  // if

        block_end = i
      }  else {
        if ( in_block && ( block_end - block_start + 1 >= 5 ) )
          blocks += new Window( block_start+1, block_end+1, consensus.substring( block_start, block_end+1 ) )
        in_block = false
        block_start = -1
        block_end = -1
      }  // if
    }  // for

    // println( "Core Blocks:" )
    // blocks foreach { case(block) =>
    //   println( block.start + "-" + block.end + " " + block.seq )
    // }  // foreach
  }  // makeBlocks

// *****************************************************************************

}  // CoreBlocks
