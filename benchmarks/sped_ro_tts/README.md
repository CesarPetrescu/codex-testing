# SpeD Romanian TTS benchmark

Controlled synthetic-speech evaluation of the Romanian SpeD Parakeet 110M checkpoint through its portable CTC-greedy ONNX export.

## Voices and services

- Microsoft Edge `ro-RO-AlinaNeural`
- Microsoft Edge `ro-RO-EmilNeural`
- Google Translate TTS via gTTS (`ro`)
- Kiprio public Romanian demo API (`ro`)
- Piper `ro_RO-mihai-medium`
- eSpeak NG Romanian baseline

## Test categories

- clean Romanian prose
- diacritics and difficult phonemes
- numbers, dates, and time
- technical infrastructure vocabulary
- Romanian-English developer code-switching
- spoken self-correction

## Output

The workflow exports:

- original and normalized WAV audio
- exact raw SpeD transcripts
- diacritic-sensitive WER/CER
- accent-insensitive WER/CER
- per-clip CPU latency, RTF, and RTFx
- CSV, JSON, Markdown, and HTML reports
- model checksums and runtime environment

## Scope limitation

This is a pronunciation and vocabulary stress test. Synthetic voices have cleaner acoustics and more regular pacing than real speakers, so these results must not be treated as production dictation accuracy. A proper decision still requires Romanian microphone recordings from multiple speakers, rooms, devices, accents, speaking rates, and spontaneous Romanian-English technical speech.

The ONNX export uses the checkpoint's CTC head with greedy decoding. It is intentionally reported separately from the original model's TDT and CTC beam-search/KenLM results.
