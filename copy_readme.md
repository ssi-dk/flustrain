# Inholdsfortegnelse

### report.txt
Denne auto-genererede rapport opsummerer alle konsensus-sekvenser, og indeholder blandt andet:
* Dybde og dækning (_coverage_)
* Identitet til tætteste reference
* Evt. problemer med alignment, assembly, eller proteiner/peptider i segmenterne.
    
### commit.txt
Denne fil indeholder en kode, angiver versionen af workflowet.

### tmp
Denne mappe indeholder interne filer fra workflowet som kan bruges til detaljeret fejlfinding.

### consensus
Indeholder konsensus-sekvenser
* Filer, der slutter på ".fna" er nukleotidsekvenser, dem på ".faa" er peptider/proteiner.
* Filer, der starter med "consensus" er alle konsensus-sekvenser. Dem med "curated" er kun de konsensus-sekvenser, der klarede alle kvalitetskontrol-trin.

### phylogeny
Automatiske fylogenetiske træer. Kun de "curated" konsensussekvenser er automatisk inkluderede her. Se her for at finde kladebestemmelsen.
