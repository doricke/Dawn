import scala.io.Source

// package DORicke
/**  Package <code>DORicke</code> defines class 
  *  <code>InputFile</code>.
  */
/** Class <code>InputFile</code> implements a text file reader object.
  *
  * @author Darell O. Ricke, Ph.D.
  *
  * @version 1.0
  */
//  This class provides an object model for an input text file.
//
//  Author:     Darrell O. Ricke, Ph.D.  (mailto: d_ricke@yahoo.com)
//  Copyright:  Copyright (c) 2011 Ricke Informatics
//  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
//  Contact:   	Ricke Informatics, 37 Pilgrim Drive, Winchester, MA 01890
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ******************************************************************************
class InputFile( var fileName: String) {
/** Create a new instance of this class with the 
  * @param fileName The name of the text file to read.
  */

  var endOfFile: Boolean = false	// end of file flag
  var contents: String = ""	      // contents of text file
  var line: String = ""           // current line of text file
  private[this] var cursor: Iterator[String] = Source.fromFile(fileName, "UTF-8").getLines()    // cursor for reading text lines from a file
  nextLine()

  // ****************************************************************************
  def openFile(): Unit = {
/** Method <code>openFile</code> opens the text file for reading.
  */
    endOfFile = false
    cursor = Source.fromFile(fileName, "UTF-8").getLines()
  }  // openFile

  // ****************************************************************************
  def nextLine(): String = {
/** Method <code>nextLine</code> reads the next line from the text file.
  *
  * @return The next line from the text file.
  */
    line = ""
    if ( cursor.hasNext )
      line = cursor.next()
    else
      endOfFile = true
    line
  }  // nextLine    

  // ****************************************************************************
  def readContents(): String = {
/** Method readContents</code> reads the entire text file into a string.
  *
  * @return The contents of the file as a string.
  */
    contents = ""
    for (line <- Source.fromFile(fileName).getLines())
      contents = contents + line + "\n"
    contents
  }  // readContents

  // ****************************************************************************

}  // class InputFile

