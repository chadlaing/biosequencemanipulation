{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}

import Protolude
import Test.Hspec
import Test.QuickCheck
import Control.Exception (evaluate)
import Bio.Sequence.Manipulation
import Bio.Core.Sequence
import Data.Char
import Data.ByteString.Lazy.Char8 as B

main :: IO ()
main = hspec $ do
  describe "Bio.Sequence.Manipulation.compChar" $ do
    it "complements any nucleotide sequence" $ do
      compChar 'A' `shouldBe` 'T'
      compChar 'v' `shouldBe` 'B'
      compChar 'X' `shouldBe` 'X'

  describe "Bio.Sequence.Manipulation.revCompSeqData" $ do
    it "reverse complements any Bio.Core.SeqData" $ do
      revCompSeqData (SeqData "AATGCN") `shouldBe` SeqData "NGCATT"

  describe "Bio.Sequence.Manipulation.revQualData" $ do
    it "reverses any Bio.Core.QualData" $ do
      revQualData (QualData "RGB%!@%") `shouldBe` QualData "%@!%BGR"

  describe "Bio.Sequence.Manipulation.revComp'" $ do
    it "adds a reverse complemented Char to the head of a growing Data.ByteString.Lazy.Char8 raw sequence" $ do
      revComp' "CCGCRAN" 'A' `shouldBe` "TCCGCRAN"


--FastQSeq{qSeqData = SeqData "ATCG"
--        ,qSeqLabel = SeqLabel "test"
--        ,qQualData = QualData "$$$$$%"
--        }