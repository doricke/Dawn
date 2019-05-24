
import scala.collection.mutable.{ArrayBuffer, Map}

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
class PdbReader( file_name: String ) extends InputFile( file_name ) {

  val chains = Map[Char, StringBuilder]()

  val residues = scala.collection.mutable.Map[Char, scala.collection.mutable.ArrayBuffer[Residue]]()

  // **************************************************************************

  // **************************************************************************
  def parseAtom( line: String ): Atom = {
    val serial = line.substring( 6, 11 )
    val name = line.substring( 17, 20 )
    val chain_id = line.charAt( 21 )
    val res_seq = line.substring( 22, 26 )
    val icode = line.charAt( 26 )
    val x = line.substring( 30, 38 )
    val y = line.substring( 38, 46 )
    val z = line.substring( 46, 54 )

    new Atom( serial, name, chain_id, res_seq, icode, x, y, z, line )
  }  // parseAtom

  // **************************************************************************
  def parsePdb() {
    var atom_blank: Atom = new Atom( "", "", ' ', "", ' ', "", "", "", "" )
    var aa_residue: Residue = new Residue( atom_blank )

    var previous_id: String = ""
    while ( endOfFile == false ) 
    {
      if ( line.substring( 0, 4 ) == "ATOM" ) {
        val residue_id = line.substring( 17, 27 )
        val atom = parseAtom( line )

        // Check for first Atom
        if ( previous_id == "" ) { 
          aa_residue = new Residue( atom )
          previous_id = residue_id
        }  // if

        // Check for new Residue
        if ( residue_id != previous_id ) {
          if ( ( residues contains atom.chain_id ) == false )
            residues += atom.chain_id -> new scala.collection.mutable.ArrayBuffer[Residue]()

          residues( aa_residue.chain_id ).append( aa_residue )
          // aa_residue.snap
          aa_residue = new Residue( atom )
          previous_id = residue_id
        }  // if
       
        aa_residue.atoms += atom 
      }  // if

      nextLine
    }  // while

    residues( aa_residue.chain_id ).append( aa_residue )
    // aa_residue.snap
   
    setChains 
  }  // method parsePdb

  // **************************************************************************
  def setChains() {
    residues foreach {case(chain, aminos) =>
      if ( ( chains contains chain ) == false )  chains( chain ) = new StringBuilder()
      aminos foreach {case(amino) => chains( chain ).append( amino.letter ) }
    }  // foreach
  }  // setChains

  // **************************************************************************
  def snap {
    chains foreach {case(chain, seq) => println( chain + "  " + seq ) }
  }  // snap

  // **************************************************************************

}  // class PdbReader
