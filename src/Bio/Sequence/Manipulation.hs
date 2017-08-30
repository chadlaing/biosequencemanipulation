{-|
Module      : Bio.Sequence.Manipulation
Description : Functions for manipulating Bio.Core.Sequence data
Copyright   : (c) Chad R Laing, 2017
License     : GPL-3
Maintainer  : chadlaing@inoutbox.com
Stability   : beta
Portability : POSIX

This module takes a Haskell-only approach to sequence manipulation
  that is clear, and builds on the Bio.Core.Sequence classes. It
  functions for Bio.Core.SeqData, Bio.Core.QualData. The module provides a
  Sequence Type that includes Amino / Nucl options for both FastA and FastQ
  data. It is to be used as a basis for future Bio.* modules dealing with
  sequence representation.
-}

{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}

module Bio.Sequence.Manipulation where


import Protolude
import Bio.Core.Sequence (unSD, unQD, SeqData(..), QualData(..), SeqLabel(..), BioSeqQual(..), Offset(..), BioSeq(..), seqdata, seqlabel, seqid, seqlength, seqqual, toFastQ, toFasta)
import Data.Char
import Data.ByteString.Lazy.Char8 as B

-- |The Type for creating either nucleotide or amino acid based sequences
data Sequence a where
  Nucl :: FastSeq -> Sequence FastSeq
  Amino :: FastSeq -> Sequence FastSeq
deriving instance Show (Sequence a)


-- | The Type for creating FastA or FastQ sequences. FastQ includes the
-- QualData for the sequence.
data FastSeq = FastA FastASeq | FastQ FastQSeq
               deriving (Eq, Show)


-- | FastA data representation
data FastASeq =
  FastASeq{
      aSeqData :: SeqData
    , aSeqLabel :: SeqLabel
  } deriving (Eq, Show)


-- | FastQ data representation
data FastQSeq =
  FastQSeq{
      qSeqData :: SeqData
    , qSeqLabel :: SeqLabel
    , qQualData :: QualData
  } deriving (Eq, Show)


-- | Implementation of Bio.Core.BioSeq for FastA data
instance BioSeq FastASeq where
  seqid = SeqLabel . B.takeWhile (/= ' ') . unSL . aSeqLabel
  seqheader = aSeqLabel
  seqdata = aSeqData
  seqlength = Offset .  B.length . unSD . aSeqData


-- | Implementation of Bio.Core.BioSeq for FastQ data
instance BioSeq FastQSeq where
  seqid = SeqLabel . B.takeWhile (/= ' ') . unSL . qSeqLabel
  seqheader = qSeqLabel
  seqdata = qSeqData
  seqlength = Offset . B.length . unSD . qSeqData


-- | Implementation of Bio.Core.BioSeqQual for FastQ data
instance BioSeqQual FastQSeq where
  seqqual = qQualData


-- |Reverse complements any Sequence. For amino acid sequences, reverse
-- complement does not make sense, so it returns the original sequence.
-- For nucleotide sequences, both FastA and FastQ get the sequence reverse
-- complemented, while FastQ additionally has the quality string reversed.
revCompSequence :: Sequence a -> Sequence a
revCompSequence x@(Amino _) = x
revCompSequence (Nucl (FastA s)) =
  Nucl . FastA $ s{aSeqData = revCompSeqData . seqdata $ s}
revCompSequence (Nucl (FastQ s)) =
  Nucl . FastQ $ s{qSeqData = revCompSeqData $ seqdata s
                  ,qQualData = revQualData . seqqual $ s}


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

