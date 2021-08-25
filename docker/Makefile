IMAGE_NAME = lr-longbow
VERSION = 0.3.1

TAG1 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):$(VERSION)
TAG2 = quay.io/broad-long-read-pipelines/$(IMAGE_NAME):latest
TAG3 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG4 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	eval `ssh-agent`; \
	ssh-add ~/.ssh/id_rsa;\
	docker build --ssh default --squash -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) . ;\
	kill $${SSH_AGENT_PID} 

build_no_cache:
	eval `ssh-agent`; \
	ssh-add ~/.ssh/id_rsa;\
	docker build --ssh default --squash --no-cache -t $(TAG1) -t $(TAG2) -t $(TAG3) -t $(TAG4) . ;\
	kill $${SSH_AGENT_PID} 

push:
	docker push $(TAG1)
	docker push $(TAG2)
	docker push $(TAG3)
	docker push $(TAG4)
