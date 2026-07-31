# CEL archive audit

Raw archives were read and extracted without modifying `data/raw/`.

 accession          archive archive_members cel_members unique_cel_names
  GSE44076 GSE44076_RAW.tar             246         246              246
  GSE41258 GSE41258_RAW.tar             390         390              390
 extracted_cels readable_cels metadata_matches status
            246           246              246   PASS
            390           390              390   PASS

All CEL filenames were unique, readable through gzip, and matched one GEO sample identifier.
