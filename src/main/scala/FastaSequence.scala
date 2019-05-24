
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

import java.io._

// import scala.collection.mutable.Map

// ****************************************************************************
class FastaSequence extends Serializable {

  // **************************************************************************
  val annotation = collection.mutable.Map[String, String]()	// sequence annotation
  var name: String = ""		// current sequence name
  var description: String = ""	// current sequence description
  var sequence: String = ""	// current sequence 
  var mid: String = ""		// barcode

  // **************************************************************************
  def getLength: Int = sequence.length

  // **************************************************************************
  def compliment( base: Char ): Char = {
    base match {
      case 'A' => 'T'
      case 'a' => 't'
      case 'C' => 'G'
      case 'c' => 'g'
      case 'G' => 'C'
      case 'g' => 'c'
      case 'T' => 'A'
      case 't' => 'a'
      case 'N' => 'N'
      case 'n' => 'n'
      case 'R' => 'Y'
      case 'r' => 'y'
      case 'Y' => 'R'
      case 'y' => 'r'
      case '.' => '.'
      case '-' => '-'
      case _ => '?'
    }  // match
  }  // compliment

  // **************************************************************************
  def parseAnnotation() = {
    val tokens = description.split( " /" )
    tokens foreach {case(token) =>
      if ( token.size > 0 ) {
        val tuple = token.split( "=" )
        if ( tuple.size > 1 )
          annotation += tuple( 0 ).filterNot( _ == '/' ) -> tuple( 1 ).filterNot( _ == '"' )
      }  // if
    }  // foreach
  }  // parseAnnotation

  // **************************************************************************
  def reverseCompliment: String = {
    var seqReverse = new StringBuilder()
    for ( i <- sequence.length-1 to 0 by -1 )
      seqReverse += compliment( sequence.charAt( i ) )
    seqReverse.toString()
  }  // reverseCompliment

  // **************************************************************************
  def parseHeader( line: String, delim: Char )
  {
    // println( "FastaSequence.parseHeader: line = " + line )

    // Check for no header line.
    if ( line.length <= 0 )
    {
      return
    }  // if

    // Check for invalid header line.
    if ( line.charAt( 0 ) != delim )
    {
      return
    }  // if

    var index1: Int = line.indexOf( ' ' )
    var index2: Int = line.indexOf( '\t' )
    if ( ( index2 > 0 ) && ( index2 < index1 ) )  index1 = index2
    if ( ( index1 < 0 ) && ( index2 > 0 ) )       index1 = index2

    if ( index1 < 0 )
      name = line.substring( 1 ).trim()
    else
    {
      name = line.substring( 1, index1 )

      if ( index1 + 1 < line.length )
        description = line.substring( index1 + 1 ).trim()
    }  // if

    parseAnnotation()
  }  // method parseHeader

  // **************************************************************************
  def parseHeader( line: String ) { parseHeader( line, '>' ) }

  // **************************************************************************
  def toBlock( str: String ): String = {
    var block: String = ""

    if ( ( str == null ) || ( str.length < 1 ) )
    { "" }
    else {
      var strStart: Int = 0
      var strEnd: Int = 50
      while ( strEnd < str.length ) {
        block += str.substring( strStart, strEnd ) + "\n"
        strStart = strEnd
        strEnd += 50
      }  // while
      
      block += str.substring( strStart ) + "\n"
      block
    }  // if
  }  // toBlock

  // **************************************************************************
  def toBlock: String = toBlock( sequence )

  // **************************************************************************
  def to_string = { ">" + name + " " + description + "\n" + toBlock }

  // **************************************************************************
}  // class FastaSequence
