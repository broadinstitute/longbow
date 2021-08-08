---
layout: default
title: Introduction
nav_order: 1
description: "Annotation and segmentation of MAS-seq data."
permalink: /
---

# Longbow

## Introduction

Longbow is a processing toolkit for **M**ultiplexed **A**rray**S** **seq**uencing (MAS-seq) long read data. It is primarily used to annotate adapter and transcript boundaries multiplexed within long reads (termed "arrays"), segment annotated reads, filter out malformed arrays, and prepare output for use with other downstream tools. 

## Features

 * Very accurate results (Al'Khafaji et al., 2021, manuscript in prep.).
 * Works well with both PacBio HiFi (rq >= 0.99) and uncorrected reads (-1 <= rq < 0.99).
 * Flexible model structure (support for different array configurations; bulk RNA and SlideSeq support in progress).
 * Easy to install
 * Open source (BSD-3-Clause License)

## Documentation

 * Online documentation: you're already here!
 * Offline documentation: available in the `docs/` subdirectory in the repository and in the downloaded .tar.gz distribution.

## Issue tracker

Please feel free to file bug reports, feature requests, etc., in [our issue tracker](https://github.com/broadinstitute/longbow/issues). We're eager to hear users' feedback!