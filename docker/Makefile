IMAGE_NAME = lr-longbow
VERSION = 0.5.0

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	eval `ssh-agent`; \
	ssh-add ~/.ssh/id_rsa;\
	docker build --ssh default --squash -t $(TAG1) -t $(TAG2) . ;\
	kill $${SSH_AGENT_PID} 

build_no_cache:
	eval `ssh-agent`; \
	ssh-add ~/.ssh/id_rsa;\
	docker build --ssh default --squash --no-cache -t $(TAG1) -t $(TAG2) . ;\
	kill $${SSH_AGENT_PID} 

push:
	docker push $(TAG1)
	docker push $(TAG2)

