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





--FastQSeq{qSeqData = SeqData "ATCG"
--        ,qSeqLabel = SeqLabel "test"
--        ,qQualData = QualData "$$$$$%"
--        }