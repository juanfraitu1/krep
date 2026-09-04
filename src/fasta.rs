//! Minimal FASTA reader and writer.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// A single FASTA entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub header: String,
    pub seq: Vec<u8>,
}

impl Record {
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
