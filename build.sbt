name := "Dawn"

version := "0.1"

scalaVersion := "2.12.6" 

resolvers += "Typesafe Repository" at "http://repo.typesafe.com/typesafe/releases/"

// mainClass in assembly := Some("Dawn.Main")

libraryDependencies ++= Seq(
    "org.slf4j" % "slf4j-simple" % "1.6.4",
    "com.typesafe.akka" % "akka-actor_2.12" % "2.5.17"
 )

