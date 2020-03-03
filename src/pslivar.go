package main

import (
	"bufio"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/biogo/hts/fai"
	"github.com/brentp/faidx"
	"github.com/brentp/gargs/process"
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func makeCommands(fasta *faidx.Faidx, args []string, directory string) chan string {
	ch := make(chan string)
	go func() {
		L := 50000000
		defer close(ch)
		recs := make([]fai.Record, 0, 32)
		for _, idx := range fasta.Index {
			recs = append(recs, idx)
		}
		sort.Slice(recs, func(i, j int) bool {
			return recs[i].Start < recs[j].Start
		})
		for k, rec := range recs {
			for i := 0; i < rec.Length; i += L {
				// on the first, split in half to get smaller size so we can get to any possible errors sooner.
				if k == 0 && i == 0 {
					Lh := L / 2
					var region = strconv.Quote(fmt.Sprintf("%s:%d-%d", rec.Name, i+1, min(i+Lh, rec.Length)))
					ch <- fmt.Sprintf("SLIVAR_QUIET= SLIVAR_SUMMARY_FILE=%s/%d-%da.summary.tsv slivar expr %s --region %s", directory, k, i, strings.Join(args, " "), region)
					time.Sleep(400 * time.Millisecond)
					region = strconv.Quote(fmt.Sprintf("%s:%d-%d", rec.Name, i+1+Lh, min(i+L, rec.Length)))
					ch <- fmt.Sprintf("SLIVAR_QUIET=yes SLIVAR_SUMMARY_FILE=%s/%d-%db.summary.tsv slivar expr %s --region %s", directory, k, i, strings.Join(args, " "), region)
				} else {
					var region = strconv.Quote(fmt.Sprintf("%s:%d-%d", rec.Name, i+1, min(i+L, rec.Length)))
					var scmd = fmt.Sprintf("SLIVAR_QUIET=yes SLIVAR_SUMMARY_FILE=%s/%d-%d.summary.tsv slivar expr %s --region %s", directory, k, i, strings.Join(args, " "), region)
					ch <- scmd
				}
			}
		}

	}()
	return ch
}

func getChromStart(cmd string) (chrom string, start int) {
	re := regexp.MustCompile(`--region ([^:]+):(\d+)`)
	m := re.FindStringSubmatch(cmd)

	start, err := strconv.Atoi(m[2])
	if err != nil {
		log.Fatal("[pslivar] error parsing region")
	}
	return m[1], start

}

type sampleCol struct {
	sample string
	col    int
}

func writeSummaryCommands(tempDir string) {
	fields := make([]string, 0, 16)
	samples := make([]string, 0, 16)

	counts := make(map[sampleCol]int, 20)

	files, err := ioutil.ReadDir(tempDir)
	if err != nil {
		panic(err)
	}
	for fi, f := range files {
		fh, err := os.Open(filepath.Join(tempDir, f.Name()))
		if err != nil {
			log.Fatal(err)
		}
		b := bufio.NewReader(fh)
		i := -1

		for {
			i++
			line, err := b.ReadString('\n')
			if err != nil && err != io.EOF {
				log.Fatal(err)
			}
			if len(line) == 0 {
				break
			}

			toks := strings.Split(strings.TrimSpace(line), "\t")
			// deal with header
			if i == 0 {
				if len(fields) == 0 {
					fields = toks[1:len(toks)]
				} else {
					for i, f := range toks[1:len(toks)] {
						if f != fields[i] {
							log.Fatalf("error fields not matching between files %s vs %s", f, fields[i])
						}
					}
				}
				continue
			}
			var sample = toks[0]

			if fi == 0 {
				samples = append(samples, sample)
			} else {
				if samples[i-1] != sample {
					log.Fatalf("samples did not match %s vs %s", samples[i-1], sample)
				}
			}

			for icol := 1; icol < len(toks); icol++ {
				key := sampleCol{sample: sample, col: icol - 1}
				val, err := strconv.Atoi(toks[icol])
				if err != nil {
					log.Fatalf("error parsing field as int: %s", err)
				}
				if n, ok := counts[key]; !ok {
					counts[key] = val
				} else {
					counts[key] = n + val
				}

			}

			if err == io.EOF {
				break
			}
		}
		fh.Close()
		os.Remove(filepath.Join(tempDir, fh.Name()))
	}

	path := os.Getenv("SLIVAR_SUMMARY_FILE")
	if path == "" {
		path = "slivar_summary.tsv"
	}
	fmt.Fprintf(os.Stderr, "writing final summary file to %s\n", path)

	out, err := os.Create(path)
	if err != nil {
		log.Fatal(err)
	}

	out.WriteString("sample\t" + strings.Join(fields, "\t"))
	out.WriteString("\n")
	for _, sample := range samples {
		cols := make([]string, 0, len(fields))
		cols = append(cols, sample)
		for ifield := range fields {
			key := sampleCol{sample: sample, col: ifield}
			cols = append(cols, strconv.Itoa(counts[key]))
		}
		out.WriteString(strings.Join(cols, "\t"))
		out.WriteString("\n")
	}
	out.Close()

}

func main() {
	args := make([]string, 0, 20)
	processes := runtime.GOMAXPROCS(0)
	var fasta *faidx.Faidx

	skipNext := false
	for i, arg := range os.Args {
		if i == 0 {
			continue
		}
		if skipNext {
			skipNext = false
			continue
		}
		if arg == "--processes" {
			var err error
			processes, err = strconv.Atoi(os.Args[i+1])
			if err != nil {
				log.Println("error parsing number of processes")
				log.Fatal(err)
			}
			skipNext = true
			continue
		}
		if arg == "-o" || arg == "--out-vcf" {
			log.Printf("[pslivar] ignoring '%s' argument as pslivar always sends VCF output to stdout", arg)
			skipNext = true
			continue
		}
		if arg == "--fasta" {
			var err error
			fasta, err = faidx.New(os.Args[i+1])
			if err != nil {
				log.Println("error parsing fasta file")
				log.Fatal(err)
			}
			skipNext = true
			continue

		}
		if arg == "expr" && i == 1 {
			continue
		}
		if arg == "--region" {
			log.Fatal("[pslivar] error can't specify --region to pslivar")
		}
		L := len(args)
		if arg[0] != '-' && L > 0 && args[L-1][0] == '-' && args[L-1] != "-" {
			arg = strconv.Quote(arg)
		}
		args = append(args, arg)
	}
	if fasta == nil {
		log.Fatal("[pslivar] run this command with the same arguments as slivar expr, but --fasta $reference is required")
	}
	runtime.GOMAXPROCS(processes)
	td, err := ioutil.TempDir("", "pslivar")
	if err != nil {
		panic("cant create tmpdir")
	}
	defer os.RemoveAll(td)
	commands := makeCommands(fasta, args, td)
	cancel := make(chan bool)
	opts := &process.Options{Ordered: true}
	out := bufio.NewWriter(os.Stdout)
	defer out.Flush()

	i := -1
	for cmd := range process.Runner(commands, cancel, opts) {
		i++
		if cmd.ExitCode() != 0 {
			if i == 0 {
				os.Stderr.Close()
			}
			log.Fatalf("[pslivar] error running command: %s", cmd.CmdStr)
			cancel <- true
		}
		chrom, start := getChromStart(cmd.CmdStr)
		if i == 0 {
			n, err := io.Copy(out, cmd)
			if err != nil {
				log.Printf("error at %s:%d after %d bytes", chrom, start, n)
				panic(err)
			}
			cmd.Close()
			out.Flush()
			continue
		}
		var done bool
		// skip the header and any variants that overlapped the previous region.
		for {
			line, err := cmd.ReadString('\n')
			if err != nil && err != io.EOF {
				panic(err)

			}
			if err == io.EOF {
				done = true
			}
			if len(line) == 0 {
				done = true
				break
			}
			if line[0] == '#' {
				continue
			}
			toks := strings.SplitN(line, "\t", 5)
			vstart, err := strconv.Atoi(toks[1])
			if err != nil {
				log.Fatal(err)
			}
			if vstart < start {
				continue
			}
			out.WriteString(line)
			break
		}
		if !done {
			_, err := io.Copy(out, cmd)
			if err != nil {
				panic(err)
			}
		}
		cmd.Close()

	}

	writeSummaryCommands(td)

}
