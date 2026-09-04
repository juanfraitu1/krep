//! Minimal FASTA reader and writer.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// The sequence identifier of a FASTA header: everything up to the first
/// whitespace. FASTA headers carry a free-text description after the ID
/// ("NC_060925.1 Homo sapiens isolate CHM13 chromosome 1, ..."), which must
/// never reach a BED/GTF column — those are whitespace-delimited formats.
pub fn seq_id(header: &str) -> &str {
    header.split_whitespace().next().unwrap_or(header)
}

/// A single FASTA entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub header: String,
    pub seq: Vec<u8>,
}

impl Record {
    /// Identifier portion of this record's header.
    pub fn id(&self) -> &str {
        seq_id(&self.header)
    }

    pub fn new(header: impl Into<String>, seq: Vec<u8>) -> Self {
        Self {
            header: header.into(),
            seq,
        }
    }
}

/// Read all records from a FASTA file, preserving case.
pub fn read<P: AsRef<Path>>(path: P) -> io::Result<Vec<Record>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(prev) = current_header.take() {
                records.push(Record::new(prev, current_seq));
                current_seq = Vec::new();
            }
            current_header = Some(header.to_string());
        } else {
            // Trim whitespace. Case is preserved so callers can distinguish
            // soft-masked (lowercase) from unmasked (uppercase) bases.
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            current_seq.extend_from_slice(trimmed.as_bytes());
        }
    }

    if let Some(header) = current_header {
        records.push(Record::new(header, current_seq));
    }

    Ok(records)
}

/// Write records to a FASTA file, wrapping sequence lines at 80 characters.
pub fn write<P: AsRef<Path>>(records: &[Record], path: P) -> io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for rec in records {
        writeln!(writer, ">{}", rec.header)?;
        for chunk in rec.seq.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Streaming FASTA reader: yields one record at a time so that a
/// multi-gigabyte genome never has to be resident in full.
///
/// `fasta::read` loads every record at once, which is fine for a single
/// chromosome but impossible for a 3.1 Gb genome on a small machine.
pub struct FastaStream<R: BufRead> {
    reader: R,
    pending_header: Option<String>,
    line: String,
    done: bool,
}

impl FastaStream<BufReader<File>> {
    /// Open a FASTA file for streaming.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(BufReader::with_capacity(1 << 20, file)))
    }
}

impl<R: BufRead> FastaStream<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            pending_header: None,
            line: String::new(),
            done: false,
        }
    }

    /// Read the next record, or `None` at end of file.
    ///
    /// `uppercase` converts the sequence in place; k-mer counting needs this
    /// so that soft-masked (lowercase) input does not change k-mer keys.
    pub fn next_record(&mut self, uppercase: bool) -> io::Result<Option<Record>> {
        if self.done {
            return Ok(None);
        }

        // Find the header for the record we are about to read.
        let header = match self.pending_header.take() {
            Some(h) => h,
            None => loop {
                self.line.clear();
                if self.reader.read_line(&mut self.line)? == 0 {
                    self.done = true;
                    return Ok(None);
                }
                if let Some(h) = self.line.trim_end().strip_prefix('>') {
                    break h.to_string();
                }
            },
        };

        let mut seq = Vec::new();
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                self.done = true;
                break;
            }
            let trimmed = self.line.trim_end();
            if let Some(h) = trimmed.strip_prefix('>') {
                self.pending_header = Some(h.to_string());
                break;
            }
            if trimmed.is_empty() {
                continue;
            }
            let start = seq.len();
            seq.extend_from_slice(trimmed.as_bytes());
            if uppercase {
                for b in &mut seq[start..] {
                    b.make_ascii_uppercase();
                }
            }
        }

        Ok(Some(Record::new(header, seq)))
    }
}

/// Streaming FASTA writer, so masked output can be emitted record by record
/// without holding the whole genome.
pub struct FastaWriter {
    writer: BufWriter<File>,
    line_width: usize,
}

impl FastaWriter {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            writer: BufWriter::with_capacity(1 << 20, File::create(path)?),
            line_width: 80,
        })
    }

    pub fn write_record(&mut self, header: &str, seq: &[u8]) -> io::Result<()> {
        writeln!(self.writer, ">{}", header)?;
        for chunk in seq.chunks(self.line_width) {
            self.writer.write_all(chunk)?;
            self.writer.write_all(b"\n")?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// Read a FASTA and convert all sequence bases to uppercase.
/// This is what k-mer counting needs, because lowercase soft masking must
/// not change the k-mer keys.
pub fn read_uppercase<P: AsRef<Path>>(path: P) -> io::Result<Vec<Record>> {
    let mut records = read(path)?;
    for rec in &mut records {
        for b in &mut rec.seq {
            *b = b.to_ascii_uppercase();
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn seq_id_strips_description() {
        assert_eq!(
            seq_id("NC_060925.1 Homo sapiens isolate CHM13 chromosome 1, alternate assembly"),
            "NC_060925.1"
        );
        assert_eq!(seq_id("chr1"), "chr1");
        assert_eq!(seq_id(""), "");
    }

    #[test]
    fn read_write_roundtrip() {
        let dir = std::env::temp_dir();
        let path = dir.join("repeatmasker-lite-test.fa");
        {
            let mut f = File::create(&path).unwrap();
            f.write_all(b">chr1\nACGT\nACGT\n>chr2\nTGCA\n").unwrap();
        }
        let recs = read(&path).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].header, "chr1");
        assert_eq!(recs[0].seq, b"ACGTACGT");
        assert_eq!(recs[1].seq, b"TGCA");

        let path2 = dir.join("repeatmasker-lite-test2.fa");
        write(&recs, &path2).unwrap();
        let recs2 = read(&path2).unwrap();
        assert_eq!(recs, recs2);
    }
}
