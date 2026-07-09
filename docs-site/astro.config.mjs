import { defineConfig } from 'astro/config';
import { unified } from '@astrojs/markdown-remark';
import starlight from '@astrojs/starlight';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';

export default defineConfig({
  site: 'https://nkzono99.github.io',
  base: '/BEACH',
  outDir: '../build/starlight-site',
  build: {
    format: 'file',
  },
  markdown: {
    processor: unified({
      remarkPlugins: [remarkMath],
      rehypePlugins: [[rehypeKatex, { strict: false, throwOnError: false }]],
    }),
  },
  integrations: [
    starlight({
      title: 'BEACH',
      defaultLocale: 'root',
      locales: {
        root: {
          label: '日本語',
          lang: 'ja',
        },
        en: {
          label: 'English',
        },
      },
      customCss: ['./src/styles/custom.css'],
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/Nkzono99/BEACH',
        },
      ],
      sidebar: [
        {
          label: 'はじめに',
          translations: { en: 'Start' },
          items: [
            { slug: 'index' },
            { slug: 'output-guide' },
          ],
        },
        {
          label: 'ユーザーガイド',
          translations: { en: 'User Guide' },
          items: [
            { slug: 'parameters' },
            { slug: 'configuration' },
            { slug: 'python-postprocess-api' },
            { slug: 'postprocess-tutorial' },
            { slug: 'configuration-recipes' },
          ],
        },
        {
          label: '数値アルゴリズム',
          translations: { en: 'Numerics' },
          items: [
            { slug: 'algorithms' },
            { slug: 'field-solvers' },
            { slug: 'particle-charge-loop' },
            { slug: 'fmm-core' },
            { slug: 'batch-duration-stability' },
          ],
        },
        {
          label: '開発者向け',
          translations: { en: 'Developers' },
          items: [
            { slug: 'workflow' },
            { slug: 'fortran-dependency-map' },
            { slug: 'agent-user-guide' },
            { label: 'Fortran API', link: 'https://nkzono99.github.io/BEACH/fortran/' },
            {
              label: 'GitHub Repository',
              link: 'https://github.com/Nkzono99/BEACH',
              attrs: { target: '_blank', rel: 'noreferrer' },
            },
          ],
        },
      ],
    }),
  ],
});
