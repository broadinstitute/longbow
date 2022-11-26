---
layout: default
title: Introduction
nav_order: 1
description: "Annotation and segmentation of MAS-seq data."
permalink: /
---

# Longbow

## Introduction

Longbow is a processing toolkit for **M**ultiplexed **A**rray**S** **seq**uencing (MAS-seq) long read data. It is primarily used to annotate adapter and transcript boundaries within multiplexed long reads (termed "arrays"), segment annotated reads, filter out malformed arrays and cDNAs, and prepare output for use with other downstream tools. 

## Features

 * Works with both low error-rate reads (PacBio HiFi `rq >= 0.99` reads) and high error-rate reads (CLR/uncorrected `-1 <= rq < 0.99` reads).
 * Flexible hierarchical model structure (can mix-and-match array and cDNA models).
 * Built-in support for several different array orderings (10-element, 15-element, and 16-element MAS-seq arrays; Iso-Seq libraries, etc.) and cDNA designs (e.g. 10x 5' single-cell kit, 10x 3' single-cell kit, bulk RNA, SlideSeq)).
 * Strong model-based filtering for arrays and cDNAs that do not conform to the expected design. 
 * Open source (BSD-3-Clause License).

## Issues

Please feel free to file bug reports, feature requests, etc., at [our issue tracker](https://github.com/broadinstitute/longbow/issues). We're always eager to hear our users' feedback!
