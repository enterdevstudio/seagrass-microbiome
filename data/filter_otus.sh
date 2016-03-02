filter_otus_from_otu_table.py -i otu_table_json_even100.biom -o otu_table_even100_host_assc.biom --negate_ids_to_exclude -e host_associated_otus.txt

biom convert -i otu_table_even100_host_assc.biom -o otu_table_json_even100_host_assc.biom  --table-type="OTU table" --to-json --header-key taxonomy




filter_otus_from_otu_table.py -i otu_table_json_even100.biom -o otu_table_even100_env_assc.biom --negate_ids_to_exclude -e env_associated_otus.txt

biom convert -i otu_table_even100_env_assc.biom -o otu_table_json_even100_env_assc.biom  --table-type="OTU table" --to-json --header-key taxonomy