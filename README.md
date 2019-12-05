# talks
NeuroData Talks

Files put in this bucket are hosted at the domain `neurodata.io/talks/*`

When a commit it made to the master branch, it triggeres an AWS CodeBuild through CodePipeline, which updates `nd-talks` S3 Bucket.

The bucket is then hosted from Cloudfront through the S3 static bucket endoint.
