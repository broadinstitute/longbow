---
layout: default
title: Introduction
nav_order: 1
description: "Annotation and segmentation of MAS-seq data."
permalink: /
---

# Longbow

## Introduction

Longbow is a processing toolkit for **M**ultiplexed **A**rray**S** **seq**uencing (MAS-seq) long read data. It is primarily used to annotate adapter and transcript boundaries within multiplexed long reads (termed "arrays"), segment annotated reads, filter out malformed arrays, and prepare output for use with other downstream tools. 

## Features

 * Works with both low error-rate reads (PacBio HiFi `rq >= 0.99` reads) and high error-rate reads (CLR/uncorrected `-1 <= rq < 0.99` reads).
 * Flexible model structure (currently supports different array orderings; bulk RNA and SlideSeq support in progress).
 * Open source (BSD-3-Clause License).

## Issues

Please feel free to file bug reports, feature requests, etc., at [our issue tracker](https://github.com/broadinstitute/longbow/issues). We're always eager to hear our users' feedback!
