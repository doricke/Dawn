
import scala.collection.mutable.{ Map }

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

class Aligner( fileName: String, target_mers: Map[String, Protein], wordSize: Int, reader: Map[String, PdbReader] ) {

  val conserved = new Conserved()

  val queries = Map[String, Protein]()
  val gene = Map[String, Map[String, Int]]()
  val family = Map[String, Map[String, Int]]()

  val geneMsa = Map[String, Msa]()
  val familyMsa = Map[String, Msa]()

  // ***************************************************************************
  def writeConservation( fasta: FastaSequence, gene_conservation: Map[Int, Float], family_conservation: Map[Int, Float], template: String, observed: Map[Int, Map[Char, Int]] ) {
    System.out.println( "Aligner writeConservation OutputFile name: " + fileName + ".csv" )

    val out = new OutputFile( fileName + ".csv" )
    val seq = fasta.sequence
    for ( i <- 1 until seq.size+1 ) {
      // val cons = if ( ( family_conservation( i ) > gene_conservation( i ) ) && ( gene_conservation( i ) > 0 ) ) family_conservation( i ) else gene_conservation( i )
      val cons = if ( template.charAt( i-1 ) != ' ' ) family_conservation( i ) else gene_conservation( i )
      out.write( i + "\t" + seq( i-1 ) + "\t" + cons + "\t" )
      val counts = observed( i ).toList.sortBy{_._2}
      val nc = conserved.nonconserved( seq(i-1), observed( i ) )
      val nc_percent = if ( counts.size > 1 ) (nc.toDouble*100.0)/counts( counts.size-1 )._2.toDouble else 0.0
      out.write( 100.0-nc_percent + "\t" + nc_percent + "\t" )
      for ( j <- counts.size-1 to 0 by -1 ) 
        out.write( counts( j )._1 + "\t" + counts( j )._2 + "\t" )
      out.write( "\n" )
    }  // for
  }  // writeConservation

  // ***************************************************************************
  def compare( fasta: FastaSequence ) {
    // println( "****** Aligner.compare called" )
    val query = new Protein( fasta, wordSize )

    geneMsa += ( fasta.name -> new Msa( fileName, fasta.name, wordSize ) )
    familyMsa += ( fasta.name -> new Msa( fileName, fasta.name, wordSize ) )
    geneMsa( fasta.name ).proteins += (fasta.name -> query)
    familyMsa( fasta.name ).proteins += (fasta.name -> query)

    gene += (fasta.name -> Map[String, Int]())
    family += (fasta.name -> Map[String, Int]())

    target_mers foreach {case (name, protein) =>
      val common_mers = query.mers & protein.mers
      val count = common_mers.size
 
      val size = if ( query.mers.size < protein.mers.size ) query.mers.size else protein.mers.size
      val percent = (count * 100) / size
      val slot = (percent / 10) * 10
      // println( fasta.name + "\t" + name + "\t" + percent + "\t" + count + "\t" + common_mers.mkString( "|" ) )

      if ( percent >= 65 ) {
        gene( fasta.name ) += (name -> count)
        geneMsa( fasta.name ).proteins += (name -> protein)
        geneMsa( fasta.name ).common_mers += (name -> common_mers)
      }
      // else
        if ( count >= 5 ) {
          family( fasta.name ) += (name -> count)
          familyMsa( fasta.name ).proteins += (name -> protein)
          familyMsa( fasta.name ).common_mers += (name -> common_mers)

          // println( "*** familyMsa: " + name )
        }  // if
    }  // foreach

    // println; println
    // println( fasta.name + "\t" + gene( fasta.name ).mkString( "|" ) )

    // println; println
    // println( fasta.name + "\t" + family( fasta.name ).mkString( "|" ) )

    // println
    // println( "Gene: " + fasta.name )
    geneMsa( fasta.name ).tallyMers()
    geneMsa( fasta.name ).bestMers()
    geneMsa( fasta.name ).setAnchors()
    geneMsa( fasta.name ).optimizeMsa()
    geneMsa( fasta.name ).writeMsa( fasta.name + "_gene" )
    // writeConservation( fasta, geneMsa( fasta.name ).conservation, geneMsa( fasta.name ).conservation, geneMsa( fasta.name ).template.toString, geneMsa( fasta.name ).observed.head._2 )

    println
    println( "Gene Family: " + fasta.name )
    familyMsa( fasta.name ).tallyMers()
    familyMsa( fasta.name ).bestMers()
    // familyMsa( fasta.name ).merCounts()
    familyMsa( fasta.name ).setAnchors()
    familyMsa( fasta.name ).optimizeMsa()
    familyMsa( fasta.name ).writeMsa( fasta.name + "_family" )
    // familyMsa( fasta.name ).sliceClique( reader, 20 )

    // println( "******* writeConservation called" )
    writeConservation( fasta, geneMsa( fasta.name ).conservation, familyMsa( fasta.name ).conservation, familyMsa( fasta.name ).template.toString, familyMsa( fasta.name ).observed.head._2 )
  }  // method compare
  
  // ***************************************************************************
  def process( fasta: FastaSequence ) = {
    compare( fasta ) 
  }  // receive

}  // class Aligner
