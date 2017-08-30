{-|
Module      : Bio.Sequence.Manipulation
Description : Functions for manipulating Bio.Core.Sequence data
Copyright   : (c) Chad R Laing, 2017
License     : GPL-3
Maintainer  : chadlaing@inoutbox.com
Stability   : stable
Portability : POSIX

This module takes a Haskell-only approach to sequence manipulation
  that is clear, and builds on the Bio.Core.Sequence classes. It
  functions for Bio.Core.SeqData, Bio.Core.QualData, and Sequence as
  implemented in <https://github.com/BioHaskell/biofasta Bio.Sequence.Fasta>
-}

{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}

module Bio.Sequence.Manipulation (
  -- * Functions
  --revCompSequence,
  revCompSeqData,
  revQualData
  ) where


import Protolude
import Bio.Core.Sequence (unSD, unQD, SeqData(..), QualData(..), SeqLabel(..), Offset(..), BioSeq(..), seqdata, seqlabel, seqid, seqlength)
import Data.Char
import Data.ByteString.Lazy.Char8 as B


data Sequence a where
  Nucl :: Fast a -> Sequence Fast
  Amino :: Fast a -> Sequence Fast
deriving instance Show (Sequence a)


data Fast a where
  FastA :: FastaSeq -> Fast FastaSeq
  FastQ :: FastqSeq -> Fast FastqSeq
deriving instance Show (Fast a)


data FastaSeq =
  FastaSeq{
      aSeqData :: SeqData
    , aSeqLabel :: SeqLabel
  } deriving (Eq, Show)


data FastqSeq =
  FastqSeq{
      qSeqData :: SeqData
    , qSeqLabel :: SeqLabel
    , qQualData :: QualData
  } deriving (Eq, Show)


defaultFasta :: Sequence Fast
defaultFasta = Nucl . FastA $ FastaSeq
                                {aSeqData = SeqData "ATCG"
                                ,aSeqLabel = SeqLabel "test"
                                }



instance BioSeq FastaSeq where
  seqid = SeqLabel . B.takeWhile (/= ' ') . unSL . aSeqLabel
  seqheader = aSeqLabel
  seqdata = aSeqData
  seqlength = Offset .  B.length . unSD . aSeqData


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
--revCompSequence :: Sequence -> Sequence
--revCompSequence (Seq l s q) = Seq l newSeq newQual
--  where
--    newSeq = revCompSeqData s
--    newQual = case q of
--      Nothing -> Nothing
--      Just x -> Just $ revQualData x


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

