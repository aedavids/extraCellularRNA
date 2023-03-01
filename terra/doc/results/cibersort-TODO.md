# slack conversation 

```
Newbie: Will WDL scatter/gather run in parallel?
Is this the best place to ask questions?
I have been using WDL for a while on terra. Works great. I love it.
I now find myself needing to run Cromwell/wdl in my universities data center.
I have a docker container that takes 3.5 days to run. To reduce the run time, I want to split
a large matrix into multiple parts, process each part in parallel then combine the results back
into a single file. The actual processing is done by a docker image/container

This seems like a good use case for scatter/gather. My question is if the parts are going to
run in parallel or not? E.G. if there are 8 parts, I want 8 docker container to start and run
in parallel. Previously when I ran Cromwell in the data center it seemed like “runTimeCPU = 8” 
in my wdl file was ignored.

Do I need the cluster admins to change some sort of configurations?

Kind regards
Andy
```

patrick Magee
```
For WDL related questions that do not relate to the specifics of cromwell execution, we have a 
dedicated and official community on slack for that. Please come and join us

As for this question in particular It really depends on the backend you are using. Typically 
scatter / gather calls (and really any call that can happen in parallel) will be submitted at 
the same time. I think there are cases where this is not true
```

aed
```
Thanks Patrick can you tell me more about backends? I planed on running the cromwell jar file 
from the cli. Is cromwell the backend?
```

patrick Magee
```
Backends are the underlying architecture that is used for running the actual tasks. Docker is 
one of those, but so is Google Life sciences, TES, HPC and others. you can check out their 
documentation to see the list of supported backends.
If your goal is only to use the CLI approach I would recommend MiniWDL. It is a lot more 
lightweight and probably fits your needs a bit better

heck it out here
chanzuckerberg/miniwdl
Workflow Description Language developer tools & local runner
Stars
134
Language
Python
```

Adam Nichols
```
I think of Cromwell more as an SDK than a product. MiniWDL is your command-line product.
```
