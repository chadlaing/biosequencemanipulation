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

  describe "Bio.Sequence.Manipulation.addCompChar" $ do
    it "adds a reverse complemented Char to the head of a growing Data.ByteString.Lazy.Char8 raw sequence" $ do
      addCompChar "CCGCRAN" 'A' `shouldBe` "TCCGCRAN"

  describe "Bio.Sequence.Manipulation.revCompSequence" $ do
    it "reverse complements a Nucl FastQ Sequence" $ do
      revCompSequence testFastQ `shouldBe` testFastQReversed
    it "reverse complements a Nucl FastA Sequence" $ do
      revCompSequence testFastA `shouldBe` testFastAReversed
    it "does not reverse complement an Amino Sequence" $ do
      revCompSequence testFastQAmino `shouldBe` testFastQAmino

  describe "Bio.Sequence.Manipulation.revCompFastASeq" $ do
    it "reverse complements a FastASeq" $ do
      revCompFastASeq testFastASeq `shouldBe` testFastASeqReversed


--Data for tests
testFastA :: Sequence FastSeq
testFastA = Nucl $ FastA $ MkFastASeq{aSeqData = SeqData "GGGTCNN"
                                     ,aSeqLabel = SeqLabel "test"}


testFastAReversed :: Sequence FastSeq
testFastAReversed  = Nucl $ FastA $ MkFastASeq{aSeqData = SeqData "NNGACCC"
                                    ,aSeqLabel = SeqLabel "test"}



testFastQ :: Sequence FastSeq
testFastQ = Nucl $ FastQ $ MkFastQSeq{qSeqData = SeqData "ATCG"
                    ,qSeqLabel = SeqLabel "test"
                    ,qQualData = QualData "$$$$$%"
                    }


testFastQAmino :: Sequence FastSeq
testFastQAmino = Amino $ FastQ $ MkFastQSeq{qSeqData = SeqData "ATCG"
                    ,qSeqLabel = SeqLabel "test"
                    ,qQualData = QualData "$$$$$%"
                    }

testFastASeq :: FastASeq
testFastASeq = MkFastASeq{aSeqData = SeqData "GGGTCNN"
                         ,aSeqLabel = SeqLabel "test"}



testFastASeqReversed :: FastASeq
testFastASeqReversed = MkFastASeq{aSeqData = SeqData "NNGACCC"
                         ,aSeqLabel = SeqLabel "test"}


testFastQReversed = Nucl $ FastQ $ MkFastQSeq{qSeqData = SeqData "CGAT"
                    ,qSeqLabel = SeqLabel "test"
                    ,qQualData = QualData "%$$$$$"
                    }

