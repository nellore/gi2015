#!/usr/bin/env bash
rail-rna align elastic -m sra_batch_13_sample_size_500.txt -i s3://rail-eu-west-1/sra_batch_13_sample_size_500_prep -o s3://rail-eu-west-1/sra_batch_13_sample_size_500_itn -a hg19 --region eu-west-1 -c 60 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price 0.60 --master-instance-bid-price 0.60 --no-consistent-view --deliverables itn --max-task-attempts 6 --ec2-key-name raileuw1
