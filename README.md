# Generating antibody escape and rbd affinity scores for Nextclade virus_properties.json

To generate escape score JSONs, run `ab_escape_json_generator.py`

These then need to be transferred to the `virus_properties.json` file in the `nextclade_data` repo

RBD binding affinity scores are generated by `rbd_json_generator.py`

Original authorship of the generator scripts isn't entirely clear, probably all of Jesse Bloom, Richard Neher and Cornelius Roemer contributed parts.