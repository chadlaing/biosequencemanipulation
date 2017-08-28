{-|
Module      : Bio.Sequence.ReverseComplement
Description : Functions for reverse complementing Bio.Core.Sequence data
                as implemented in Bio.Sequence.Fasta (https://github.com/BioHaskell/biofasta)
Copyright   : (c) Chad R Laing, 2017
License     : GPL-3
Maintainer  : chadlaing@inoutbox.com
Stability   : stable
Portability : POSIX

This module takes a Haskell-only approach to reverse complement a sequence
  that is clear, rather than an approach from the Benchmark Game. It provides
  functions for Bio.Core.SeqData, Bio.Core.QualData, and Sequence as
  implemented in <https://github.com/BioHaskell/biofasta Bio.Sequence.Fasta>
-}

{-# LANGUAGE NoImplicitPrelude #-}

module Bio.Sequence.ReverseComplement (
  -- * Functions
  revCompSequence,
  revCompSeqData,
  revQualData
  ) where

import Protolude
import Bio.Core.Sequence (unSD, unQD, SeqData(..), QualData(..))
import Bio.Sequence.Fasta (readFasta, toStr, Sequence(..))
import Data.Char
import Data.ByteString.Lazy.Char8 as B

-- |Complements the nucleotide base, including ambiguous characters. This
-- function is not exported.
compChar :: Char -> Char
compChar c = case c of
  'A' -> 'T'
  'a' -> 'T'
  'C' -> 'G'
  'c' -> 'G'
  'G' -> 'C'
  'g' -> 'C'
  'T' -> 'A'
  't' -> 'A'
  'U' -> 'A'
  'u' -> 'A'
  'M' -> 'K'
  'm' -> 'K'
  'R' -> 'Y'
  'r' -> 'Y'
  'Y' -> 'R'
  'y' -> 'R'
  'K' -> 'M'
  'k' -> 'M'
  'V' -> 'B'
  'v' -> 'B'
  'H' -> 'D'
  'h' -> 'D'
  'D' -> 'H'
  'd' -> 'H'
  'B' -> 'V'
  'b' -> 'V'
  x -> x

-- | The reverse complement function for 'Bio.Fasta.Sequence' data. As quality
-- is an optional component of 'Bio.Fasta.Sequence', reverse the quality data
-- so that it matches the reverse complemented sequence. Uses a strict left
-- fold for the data accumulation. 'Sequence' as implemented
-- in <https://github.com/BioHaskell/biofasta Bio.Fasta.Sequence> is as follows:
--
-- > data Sequence = Seq SeqLabel SeqData (Maybe QualData)
-- >      deriving (Show, Eq)
-- For Example:
--
-- > revCompSequence $ Seq (SeqLabel "test") (SeqData "ATGCCG") (Just $ QualData "(::::@")
-- Gives:
--
-- > Seq (SeqLabel {unSL = "test"}) (SeqData {unSD = "CGGCAT"}) (Just (QualData {unQD = "@::::("}))
--
revCompSequence :: Sequence -> Sequence
revCompSequence (Seq l s q) = Seq l newSeq newQual
  where
    newSeq = revCompSeqData s
    newQual = case q of
      Nothing -> Nothing
      Just x -> Just $ revQualData x


-- | The reverse complement function for 'Bio.Core.SeqData'
--
-- For example:
--
-- > revCompSeqData $ SeqData "ATGCCG"
-- Gives:
--
-- > SeqData {unSD = "CGGCAT"}
revCompSeqData :: SeqData -> SeqData
revCompSeqData (SeqData s) = SeqData $ B.foldl' revComp' B.empty s


-- | The reverse function for 'Bio.Core.QualData'
--
-- For example:
--
-- > revQualData $ QualData "(::::@"
-- Gives:
--
-- > QualData {unQD = "@::::("}
revQualData :: QualData -> QualData
revQualData (QualData q) = QualData $ B.reverse q


-- | The reverse complement function for 'Data.Bytestring.Char8' raw sequence
-- data.
revComp' :: B.ByteString
         -> Char
         -> B.ByteString
revComp' s c = B.cons (compChar c) s